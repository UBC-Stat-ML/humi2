#!/usr/bin/env nextflow

import java.util.stream.Collectors


pwd=new File(".").getAbsolutePath()

datasets = Channel.fromPath( 'data/*', type: 'dir' ).filter{ 
  result = !it.getName().startsWith(".") 
  if (result) println("Queuing data " + it.toString())
  return result
}.map{it.getName().toString()}

params.model = "BNB LocalLambdaMixBNB MixBNB MixNB MixYS NB Poi YS MultiHit"

params.nScans = 1000
params.nInitParticles = 10 // increase this if model initialization fails (can happen in complex mixture models with vague priors)
params.nTargets = "INF" // use this to do inference on a subset of targets (e.g. for dry runs)
params.nChains = 1
params.nCores = 1
params.onlyComputeEstimates = true

deliverableDir = 'deliverables/' + workflow.scriptName.replace('.nf','') + "_" + params.nScans + "_" + params.nInitParticles + "_" + params.nTargets + "_" + params.nChains + "/"
runsDir = deliverableDir + "runs" 

models = Arrays.asList(params.model.split("\\s+")).stream().map{
  result = "humi.models." + it
  println("Queuing model " + result)
  return result
}.collect(Collectors.toList())

process buildCode {
  cache true 
  executor 'local'
  input:
    val gitRepoName from 'nowellpack'
    val gitUser from 'UBC-Stat-ML'
    val codeRevision from '14401b92c277ddc5994ee869fc336f8f091f7374'
    val snapshotPath from "${System.getProperty('user.home')}/w/nowellpack"
  output:
    file 'code' into code
  script:
    template 'buildRepo.sh' 
}

process run {

  time '30h'
  errorStrategy 'ignore'
  cpus params.nCores

  input:
    file code
    each model from models
    each dataset from datasets
  output:
    file "${dataset}_${model}" into runs
  publishDir runsDir, mode: 'link'
"""
  ln -s $pwd/data/$dataset/final.csv ${dataset}.csv  # o.w. spark crash on it
  java -cp code/lib/\\* -Xmx1g $model   \
           --experimentConfigs.resultsHTMLPage false \
           --experimentConfigs.tabularWriter.compressed true \
           --model.initialPopCounts.dataSource $pwd/data/$dataset/initial.csv \
           --model.initialPopCounts.name counts  \
           --model.data.source ${dataset}.csv \
           --model.data.genes.name gene     \
           --model.data.targets.name sgRNA     \
           --model.data.targets.maxSize $params.nTargets \
           --model.data.experiments.name dataset     \
           --model.data.experiments.maxSize 1 \
           --model.data.histograms.name histogram     \
           --engine.nScans $params.nScans   \
           --engine.nChains $params.nChains \
           --engine.nPassesPerScan 1 \
           --engine.thinning 1 \
           --engine.nThreads Max     \
           --engine.scmInit.nParticles $params.nInitParticles \
           --engine.scmInit.temperatureSchedule.threshold 0.6 \
           --engine.scmInit.nThreads Max \
           --postProcessor humi.HumiPostProcessor \
           --postProcessor.data.targets.name sgRNA \
           --postProcessor.data.genes.name gene \
           --postProcessor.data.experiments.name dataset \
           --postProcessor.data.histograms.name histogram \
           --postProcessor.runPxviz false \
           --postProcessor.onlyComputeEstimates $params.onlyComputeEstimates
  mv results/all/`ls results/all` ${dataset}_${model}
  """
}


process analysisCode {
  executor 'local'
  input:
    val gitRepoName from 'nedry'
    val gitUser from 'UBC-Stat-ML'
    val codeRevision from 'a9abcc40abcfb285588cc4c312d8ecc0bbdad06e'
    val snapshotPath from "${System.getProperty('user.home')}/w/nedry"
  output:
    file 'code' into analysisCode
  script:
    template 'buildRepo.sh'
}

process aggregate {
  input:
    file analysisCode
    file 'exec_*' from runs.toList()
  output:
    file 'output' into aggregated
  """
  code/bin/aggregate \
    --experimentConfigs.managedExecutionFolder false \
    --experimentConfigs.resultsHTMLPage false \
    --dataPathInEachExecFolder gof.csv.gz estimates.csv.gz monitoring/logNormalizationContantProgress.csv.gz \
    --keys model.data.source as data model from arguments.tsv
  mkdir output
  mv gof.csv output
  mv estimates.csv output
  mv logNormalizationContantProgress.csv output || true
  """
}

process plot {
  scratch false
  container 'cgrlab/tidyverse'
  input:
    file aggregated
  output:
    file '*.csv'
    file '*.pdf'
  publishDir deliverableDir, mode: 'copy', overwrite: true
  """
  #!/usr/bin/env Rscript
  require("ggplot2")
  require("stringr")

  data <- read.csv("$aggregated/gof.csv")
  
  data\$model <- str_replace_all(data\$model, "[\$].*", "")
  data\$model <- str_replace_all(data\$model, "humi[.]models[.]", "")
  
  write.csv(data, file="gof-data.csv")
  
  p <- ggplot(data, aes(x = data, y = actualCoverage, colour = model, group = model)) + 
    geom_line() + 
    facet_grid(gofStatistic ~ ., scales="free") + 
    geom_hline(yintercept=0.9) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45,hjust = 1)) 

  ggsave(plot = p, filename = "gof.pdf")
  
  p <- ggplot(data, aes(x = data, y = width, colour = model, group = model)) + 
    geom_line() + 
    facet_grid(gofStatistic ~ ., scales="free") + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45,hjust = 1)) 

  ggsave(plot = p, filename = "width.pdf")
  
  
  data <- read.csv("$aggregated/estimates.csv")
  
  data\$model <- str_replace_all(data\$model, "[\$].*", "")
  data\$model <- str_replace_all(data\$model, "humi[.]models[.]", "")
  
  p <- ggplot(data, aes(x = factor(model), y = logRatio, colour = model)) + 
    coord_flip() + 
    geom_errorbar(aes(ymin=logRatioLeftBound, ymax=logRatioRightBound)) +
    geom_point() + 
    facet_grid(gene + sgrna ~ data) +
    theme_bw() + 
    xlab("Gene") + 
    ylab("log(ratio)") + 
    ggtitle("Ratio of clone sizes relative to controls", subtitle = "Bayesian hierarchical model credible intervals") + 
    geom_hline(yintercept=0) + 
    theme(legend.position="none") 
  ggsave(plot = p, filename = "intervals-multi.pdf", height = 100, width = 20, limitsize = FALSE)
  
  if (file.exists("$aggregated/logNormalizationContantProgress.csv")) {
    data <- read.csv("$aggregated/logNormalizationContantProgress.csv")
  
    data\$model <- str_replace_all(data\$model, "[\$].*", "")
    data\$model <- str_replace_all(data\$model, "humi[.]models[.]", "")
  
    write.csv(data, file="evidence.csv")
  
    #data <- data %>% filter(round > 3)
  
    p <- ggplot(data, aes(x = round, y = value, colour = model, group = model)) + 
      geom_line() + 
      facet_grid(data ~ ., scales="free") + 
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45,hjust = 1)) 

    ggsave(plot = p, filename = "evidence.pdf", height = 20)
  }

  """
}


