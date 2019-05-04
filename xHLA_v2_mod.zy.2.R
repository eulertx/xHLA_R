library(jsonlite)
library(lpSolve)
library(dplyr)
library(tidyr)


message("checking installation of reference file")
if(!file.exists("chr6/chr6.fa.pac")){system("bwa index -a bwtsw chr6/chr6.fa")}

# message("creating diamond file")
# system("diamond makedb --in data/hla.faa -d data/hla")

# message("unzipping any gz files")
# system("gunzip input_test/*.gz")
a <-list.files(path = "input_test/",pattern = "fastq.gz")
a<-as.data.frame(a)
b <- separate(a,col="a",into=c('id','SinSeq','Read','others'),sep="_",remove=FALSE)
b[b$Read=='R1',] -> bR1
b[b$Read=='R2',] -> bR2
colnames(bR1)[1] <- 'readone'
colnames(bR2)[1] <- 'readtwo'
left_join(bR1%>%dplyr::select(id, readone), bR2%>%dplyr::select(id,readtwo)) -> samples

#b<-as.data.frame(regexpr(pattern = "*_R*",text = a$a))
#samples<-substr(a$a,start = 1,stop = (b[,1]))
#samples<-as.data.frame(unique(samples))
#names(samples)<-c("id")


#samples$readone<-paste(samples$id,"_R1_001.fastq.gz.trimmed.gz",sep="")
#samples$readtwo<-paste(samples$id,"_R2_001.fastq.gz.trimmed.gz",sep="")
if(length(which(samples$id==""))>0){samples<-samples[-which(samples$id==""),]}



message("running HLA typing analysis on all specimens")

for (i in 1:dim(samples)[1])
{

#  message(paste("starting analysis ",samples$readone[i],sep=""))
#  message("getting sampled data set from fastq files")
#  seeder<-sample(1:1e6, 1)

#  system(paste("seqtk sample -s", seeder ," input_test/", samples$readone[i]," 15000 > input_test/input_test_tmp1.fastq",sep=""))
#  system(paste("seqtk sample -s", seeder ," input_test/", samples$readtwo[i]," 15000 > input_test/input_test_tmp2.fastq",sep=""))

message("Running bwa to compare chr6 to fastq with speedseq mod")
# system(paste("/gpfs/bin/common/bwa mem -t 20 chr6/chr6.fa input_test/" , samples$readone[i] , " input_test/", samples$readtwo[i] , " > input_test/",samples$id[i],".sam",sep=""))

# message("Running samtools to index and align reads")

# system(paste("samtools view -bT chr6/chr6.fa input_test/",samples$id[i],".sam > input_test/",samples$id[i],".bam",sep=""))
# system(paste("samtools sort input_test/",samples$id[i],".bam > input_test/",samples$id[i],".sorted.bam",sep=""))
# system(paste("samtools index input_test/",samples$id[i],".sorted.bam input_test/",samples$id[i],".sorted.bai",sep=""))

system(paste("speedseq align -t 20 -K /gpfs/bin/scripts/somatic_pipeline_branch/somatic_pipeline/src/speedseq.config -o input_test/", samples$id[i],".sorted.bam ", " -R \"@RG:\\tID:TEST\\tSM:TEST\\tPL:ILLUMINA\\tCN:Euler\\tLB:test\\tPU:test\" chr6/chr6.fa input_test/", samples$readone[i] , " input_test/", samples$readtwo[i], sep=""))

# paste("speedseq align -t 20 -K /gpfs/bin/scripts/somatic_pipeline_branch/somatic_pipeline/src/speedseq.config -o input_test/", samples$id[i],".sorted.bam ", " -R \"@RG:\\tID:TEST\\tSM:TEST\\tPL:ILLUMINA\\tCN:Euler\\tLB:test\\tPU:test\" chr6/chr6.fa input_test/", samples$readone[i] , " input_test/", samples$readtwo[i], sep="")



 
# speedseq align -K $base/src/speedseq.config -t $nt -o $h{$sample}{$type}{'name'}.$type.$i -R \"\@RG:\\tID:$h{$sample}{$type}{'name'}.$i\\tSM:$h{$sample}{$type}{'name'}\\tPL:ILLUMINA\\tCN:Euler\\tLB:$sample.$h{$sample}{$type}{'name'}.LIB$i\\tPU:test_pu\" $ref $r1 $r2

message("running xHLA typer.sh")

system(paste("perl bin/typer.sh input_test/",samples$id[i],".sorted.bam ",samples$id[i],sep=""))


message("reading results data")  
if (length(list.files(paste("hla-",samples$id[i],"/",sep=""),pattern=".json"))>0)
{  tmp1<-fromJSON(txt=paste("hla-",samples$id[i],"/",samples$id[i],".json",sep=""))
  write_json(x = tmp1,path = paste("output/",samples$id[i],".json",sep=""))
}


  #prepare format for first specimen
 if (samples$id[i] %in% gsub(".json","",(list.files("output/",pattern=".json"))))
{
  tmp<-fromJSON(txt=paste("output/",samples$id[i],".json",sep=""))
  sampleID<-tmp$subject_id
  outputreport<-as.data.frame(sampleID,stringsAsFactors = F)
  outputreport$analysis.date<-tmp$creation_time
  outputreport$report.version<-tmp$report_version
  outputreport$report_type<-tmp$report_type

  outputreport$hla.Ai<-tmp$hla$alleles[grep(pattern = "^A",x = tmp$hla$alleles)[1]]
  outputreport$hla.Aii<-tmp$hla$alleles[grep(pattern = "^A",x = tmp$hla$alleles)[2]]
  outputreport$hla.Bi<-tmp$hla$alleles[grep(pattern = "^B",x = tmp$hla$alleles)[1]]
  outputreport$hla.Bii<-tmp$hla$alleles[grep(pattern = "^B",x = tmp$hla$alleles)[2]]
  outputreport$hla.Ci<-tmp$hla$alleles[grep(pattern = "^C",x = tmp$hla$alleles)[1]]
  outputreport$hla.Cii<-tmp$hla$alleles[grep(pattern = "^C",x = tmp$hla$alleles)[2]]
  outputreport$hla.DPB1i<-tmp$hla$alleles[grep(pattern = "^DP",x = tmp$hla$alleles)[1]]
  outputreport$hla.DPB1ii<-tmp$hla$alleles[grep(pattern = "^DP",x = tmp$hla$alleles)[2]]
  outputreport$hla.DQB1i<-tmp$hla$alleles[grep(pattern = "^DQ",x = tmp$hla$alleles)[1]]
  outputreport$hla.DQB1ii<-tmp$hla$alleles[grep(pattern = "^DQ",x = tmp$hla$alleles)[2]]
  outputreport$hla.DRB1i<-tmp$hla$alleles[grep(pattern = "^DR",x = tmp$hla$alleles)[1]]
  outputreport$hla.DRB1ii<-tmp$hla$alleles[grep(pattern = "^DR",x = tmp$hla$alleles)[2]]
  #outputreport$crosschecked<-paste(tmp$hla$matches,collapse = ",")
  rm(tmp)
  data_out<-outputreport

  message("writing results to output/000_alldata_out.txt")
  if(file.exists("output/000_alldata_out.txt")){write.table(outputreport,append = T,file = "output/000_alldata_out.txt",col.names = F,row.names = F,sep = "\t",quote=F)}
  if(!file.exists("output/000_alldata_out.txt")){write.table(outputreport,append = F,file = "output/000_alldata_out.txt",col.names = T,row.names = F,sep = "\t",quote=F)}
}
	
  system("rm input_test/input_test_tmp*.fastq")
  system("rm input_test/*.sam")
  system("rm input_test/*.bam")
  system("rm input_test/*.bai")
  system(paste ("rm -rf hla-",samples$id[i],"/",sep="")) 
  message("done")
}

# message("remove fastq files")
# system("rm input_test/*fastq")
