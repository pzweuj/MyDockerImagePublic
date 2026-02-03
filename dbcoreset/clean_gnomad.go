package main

import (
	"bufio"
	"flag"
	"fmt"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"strings"
	"sync"
)

// 保留的 INFO 字段（用户指定的具体字段）
var keepInfo = map[string]bool{
	// _joint 后缀
	"AF_joint": true, "grpmax_joint": true, "nhomalt_joint": true, "AC_joint": true, "AN_joint": true,
	// _joint_XX 后缀
	"AF_joint_XX": true, "nhomalt_joint_XX": true, "AC_joint_XX": true, "AN_joint_XX": true,
	// _joint_XY 后缀
	"AF_joint_XY": true, "nhomalt_joint_XY": true, "AC_joint_XY": true, "AN_joint_XY": true,
	// _joint_eas 后缀
	"AF_joint_eas": true,
}

// 染色体列表
var chromosomes = []string{"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"}

type Result struct {
	Chrom       string
	LineCount   int
	PassedCount int
	Error       error
}

// 串行模式标志
var serialMode bool

func main() {
	inputDir := flag.String("i", "downloads", "输入目录")
	outputDir := flag.String("o", "output", "输出目录")
	workers := flag.Int("w", runtime.NumCPU(), "并行数")
	keepTemp := flag.Bool("keep-temp", false, "保留中间文件")
	flag.BoolVar(&serialMode, "serial", false, "串行模式：逐个染色体处理并立即合并（节省磁盘空间）")
	flag.Parse()

	fmt.Println("=" + strings.Repeat("=", 59))
	fmt.Println("GnomAD v4.1 数据库清洗 (Go 版本)")
	fmt.Println("=" + strings.Repeat("=", 59))
	fmt.Printf("输入目录: %s\n", *inputDir)
	fmt.Printf("输出目录: %s\n", *outputDir)
	fmt.Printf("并行数: %d\n", *workers)
	if serialMode {
		fmt.Println("模式: 串行（节省磁盘空间）")
	} else {
		fmt.Println("模式: 并行")
	}

	// 查找输入文件
	files := findInputFiles(*inputDir)
	if len(files) == 0 {
		fmt.Printf("错误: 在 %s 中未找到 GnomAD 文件\n", *inputDir)
		os.Exit(1)
	}
	fmt.Printf("找到 %d 个染色体文件\n", len(files))

	// 创建输出目录
	os.RemoveAll(*outputDir)
	os.MkdirAll(*outputDir, 0755)
	outputFile := filepath.Join(*outputDir, "gnomad.v4.1.slim.vcf.gz")

	if serialMode {
		// 串行模式：逐个处理，处理完立即合并并删除临时文件
		runSerialMode(files, *inputDir, *outputDir, outputFile)
	} else {
		// 并行模式：先处理完所有染色体，再合并
		runParallelMode(files, *inputDir, *outputDir, outputFile, *workers, *keepTemp)
	}
}

func findInputFiles(dir string) map[string]string {
	files := make(map[string]string)
	for _, chrom := range chromosomes {
		for _, pattern := range []string{
			"gnomad.joint.v4.1.sites.chr" + chrom + ".vcf.bgz",
			"gnomad.joint.v4.1.sites.chr" + chrom + ".vcf.gz",
		} {
			path := filepath.Join(dir, pattern)
			if _, err := os.Stat(path); err == nil {
				files[chrom] = path
				break
			}
		}
	}
	return files
}

func worker(jobs <-chan string, results chan<- Result, inputDir, tempDir string) {
	for chrom := range jobs {
		res := processChromosome(chrom, inputDir, tempDir)
		results <- res
	}
}

// 串行模式：逐个染色体处理，立即合并
func runSerialMode(files map[string]string, inputDir, outputDir, outputFile string) {
	tempFile := filepath.Join(outputDir, "temp_chr.vcf.gz")
	firstFile := true
	var totalPassed int

	for _, chrom := range chromosomes {
		if _, ok := files[chrom]; !ok {
			continue
		}

		fmt.Printf("\n处理 chr%s...\n", chrom)

		// 处理单个染色体
		res := processChromosome(chrom, inputDir, outputDir)
		if res.Error != nil {
			fmt.Printf("错误: %v\n", res.Error)
			continue
		}

		chromFile := filepath.Join(outputDir, "gnomad.chr"+chrom+".vcf.gz")
		fmt.Printf("chr%s: %d PASS\n", chrom, res.PassedCount)
		totalPassed += res.PassedCount

		if firstFile {
			// 第一个文件：直接重命名
			os.Rename(chromFile, outputFile)
			os.Rename(chromFile+".tbi", outputFile+".tbi")
			firstFile = false
		} else {
			// 后续文件：追加合并到最终文件
			fmt.Printf("  追加 chr%s 到最终文件...\n", chrom)
			if err := appendVCF(outputFile, chromFile); err != nil {
				fmt.Printf("追加失败: %v，尝试备用方案...\n", err)
				// 备用方案：重新合并
				mergeVCFs(map[string]string{
					"prev": outputFile,
					chrom: chromFile,
				}, outputFile+".new")
				os.Remove(outputFile)
				os.Remove(outputFile + ".tbi")
				os.Rename(outputFile+".new", outputFile)
			}
			// 删除临时染色体文件
			os.Remove(chromFile)
			os.Remove(chromFile + ".tbi")
		}
	}

	// 重新索引最终文件
	cmdIndex := exec.Command("bcftools", "index", "-t", outputFile)
	cmdIndex.Run()

	fmt.Println("\n完成!")
	fmt.Printf("输出: %s\n", outputFile)
	fmt.Printf("索引: %s.tbi\n", outputFile)
	fmt.Printf("总变异: %d\n", totalPassed)
}

// 并行模式：先处理完所有染色体，再合并
func runParallelMode(files map[string]string, inputDir, outputDir, outputFile string, workers int, keepTemp bool) {
	tempDir := filepath.Join(outputDir, "temp")
	os.MkdirAll(tempDir, 0755)

	jobs := make(chan string, len(files))
	results := make(chan Result, len(files))

	var wg sync.WaitGroup
	for i := 0; i < workers; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			worker(jobs, results, inputDir, tempDir)
		}()
	}

	// 发送任务
	go func() {
		for chrom := range files {
			jobs <- chrom
		}
		close(jobs)
	}()

	// 收集结果
	go func() {
		wg.Wait()
		close(results)
	}()

	// 处理结果
	filteredFiles := make(map[string]string)
	var totalPassed int
	for res := range results {
		if res.Error != nil {
			fmt.Printf("错误: %v\n", res.Error)
			continue
		}
		fmt.Printf("chr%s: %d PASS\n", res.Chrom, res.PassedCount)
		filteredFiles[res.Chrom] = filepath.Join(tempDir, "gnomad.chr"+res.Chrom+".vcf.gz")
		totalPassed += res.PassedCount
	}

	// 合并文件
	fmt.Printf("\n合并... (%d 总变异)\n", totalPassed)
	mergeVCFs(filteredFiles, outputFile)

	// 清理
	if !keepTemp {
		os.RemoveAll(tempDir)
	}

	fmt.Println("\n完成!")
	fmt.Printf("输出: %s\n", outputFile)
	fmt.Printf("索引: %s.tbi\n", outputFile)
}

func processChromosome(chrom, inputDir, tempDir string) Result {
	inputFile := ""
	for _, pattern := range []string{
		"gnomad.joint.v4.1.sites.chr" + chrom + ".vcf.bgz",
		"gnomad.joint.v4.1.sites.chr" + chrom + ".vcf.gz",
	} {
		path := filepath.Join(inputDir, pattern)
		if _, err := os.Stat(path); err == nil {
			inputFile = path
			break
		}
	}

	if inputFile == "" {
		return Result{Chrom: chrom, Error: fmt.Errorf("文件未找到")}
	}

	outputFile := filepath.Join(tempDir, "gnomad.chr"+chrom+".vcf.gz")
	passedCount, err := filterVCF(inputFile, outputFile)
	if err != nil {
		return Result{Chrom: chrom, Error: err}
	}

	// 创建索引
	cmd := exec.Command("bcftools", "index", "-t", outputFile)
	if err := cmd.Run(); err != nil {
		return Result{Chrom: chrom, Error: fmt.Errorf("索引失败: %v", err)}
	}

	return Result{Chrom: chrom, PassedCount: passedCount}
}

func filterVCF(inputFile, outputFile string) (int, error) {
	// 方案: 使用 bcftools 解压，然后用 Go 过滤，直接管道输出到 bgzip 压缩
	// bcftools view | Go 过滤 | bgzip > output.gz

	// 启动 bcftools 解压
	cmdView := exec.Command("bcftools", "view", inputFile)
	stdout, err := cmdView.StdoutPipe()
	if err != nil {
		return 0, err
	}

	// 启动 bgzip 压缩（输入来自 stdin）
	cmdGzip := exec.Command("bgzip", "-c")
	gzipStdin, err := cmdGzip.StdinPipe()
	if err != nil {
		return 0, err
	}

	// 输出到临时文件，完成后重命名
	tmpFile := outputFile + ".tmp"
	out, err := os.Create(tmpFile)
	if err != nil {
		return 0, err
	}

	if err := cmdView.Start(); err != nil {
		return 0, fmt.Errorf("启动 bcftools 失败: %v", err)
	}
	if err := cmdGzip.Start(); err != nil {
		return 0, fmt.Errorf("启动 bgzip 失败: %v", err)
	}

	// Go 过滤并直接管道给 bgzip
	reader := bufio.NewReader(stdout)
	inHeader := true
	passedCount := 0
	buf := make([]byte, 64*1024) // 64KB buffer

	for {
		line, err := reader.ReadString('\n')
		if err == io.EOF {
			break
		}
		line = strings.TrimRight(line, "\r\n")

		if inHeader {
			if strings.HasPrefix(line, "#") {
				if strings.HasPrefix(line, "##contig=") && !strings.Contains(line, "ID=chr") {
					line = strings.Replace(line, "ID=", "ID=chr", 1)
				}
				gzipStdin.Write([]byte(line + "\n"))
			} else {
				inHeader = false
			}
			continue
		}

		parts := strings.Split(line, "\t")
		if len(parts) < 8 {
			continue
		}

		if parts[6] != "PASS" {
			continue
		}

		parts[7] = filterInfo(parts[7])

		if !strings.HasPrefix(parts[0], "chr") {
			parts[0] = "chr" + parts[0]
		}

		gzipStdin.Write([]byte(strings.Join(parts, "\t") + "\n"))
		passedCount++
	}

	gzipStdin.Close()
	cmdView.Wait()

	// 直接从 bgzip stdout 复制到文件
	io.Copy(out, bufio.NewReader(cmdGzip.Stdout))
	cmdGzip.Wait()
	out.Close()

	// 重命名临时文件
	os.Rename(tmpFile, outputFile)

	return passedCount, nil
}

func mergeVCFs(filteredFiles map[string]string, outputFile string) {
	// 合并策略：直接使用 bcftools concat
	var vcfPaths []string
	for _, chrom := range chromosomes {
		if path, ok := filteredFiles[chrom]; ok {
			vcfPaths = append(vcfPaths, path)
		}
	}

	if len(vcfPaths) == 1 {
		// 只有一个文件，直接重命名
		os.Rename(vcfPaths[0], outputFile)
		os.Rename(vcfPaths[0]+".tbi", outputFile+".tbi")
		return
	}

	// 使用 bcftools concat 合并
	args := append([]string{"concat", "-Oz", "-o", outputFile}, vcfPaths...)
	cmd := exec.Command("bcftools", args...)
	if err := cmd.Run(); err != nil {
		fmt.Printf("bcftools concat 失败，尝试备用方案: %v\n", err)
		mergeVCFFallback(vcfPaths, outputFile)
		return
	}

	// 创建索引
	cmdIndex := exec.Command("bcftools", "index", "-t", outputFile)
	cmdIndex.Run()
}

func mergeVCFFallback(vcfPaths []string, outputFile string) {
	// 备用合并方案：直接拼接并管道压缩
	// 适用于已 bgzip 压缩的文件

	// 启动 bgzip 压缩
	cmdGzip := exec.Command("bgzip", "-c")
	gzipStdin, _ := cmdGzip.StdinPipe()
	gzipStdout, _ := cmdGzip.StdoutPipe()

	tmpFile := outputFile + ".tmp"
	out, _ := os.Create(tmpFile)

	cmdGzip.Start()

	// 并发：写入 bgzip 和复制到文件
	done := make(chan bool, 1)
	go func() {
		io.Copy(out, bufio.NewReader(gzipStdout))
		done <- true
	}()

	// 合并数据
	for i, path := range vcfPaths {
		// 获取表头
		cmd := exec.Command("bcftools", "view", "-h", path)
		output, _ := cmd.Output()

		if i == 0 {
			gzipStdin.Write(output)
			gzipStdin.Write([]byte("\n"))
		}

		// 解压数据行并复制
		cmd = exec.Command("bcftools", "view", path)
		stdout, _ := cmd.StdoutPipe()
		cmd.Start()
		io.Copy(gzipStdin, bufio.NewReader(stdout))
		cmd.Wait()

		fmt.Printf("  合并进度: %d/%d\n", i+1, len(vcfPaths))
	}

	gzipStdin.Close()
	<-done
	cmdGzip.Wait()
	out.Close()

	os.Rename(tmpFile, outputFile)
}

func filterInfo(info string) string {
	if info == "." || info == "" {
		return "."
	}

	parts := strings.Split(info, ";")
	var keep []string

	for _, part := range parts {
		idx := strings.Index(part, "=")
		if idx > 0 {
			key := part[:idx]
			if keepInfo[key] {
				keep = append(keep, part)
			}
		}
	}

	if len(keep) == 0 {
		return "."
	}
	return strings.Join(keep, ";")
}

func formatNumber(n int) string {
	if n >= 1000000 {
		return fmt.Sprintf("%.1fM", float64(n)/1000000)
	}
	if n >= 1000 {
		return fmt.Sprintf("%.1fK", float64(n)/1000)
	}
	return fmt.Sprintf("%d", n)
}

// appendVCF 将 srcVCF 的数据行追加到 destVCF（跳过表头）
func appendVCF(destVCF, srcVCF string) error {
	// 使用 bcftools concat 追加（会去除重复表头）
	tmpFile := destVCF + ".concat.tmp"
	args := []string{"concat", "-Oz", "-o", tmpFile, destVCF, srcVCF}
	cmd := exec.Command("bcftools", args...)
	if err := cmd.Run(); err != nil {
		return err
	}
	// 替换原文件
	os.Remove(destVCF)
	os.Remove(destVCF + ".tbi")
	os.Rename(tmpFile, destVCF)
	return nil
}
