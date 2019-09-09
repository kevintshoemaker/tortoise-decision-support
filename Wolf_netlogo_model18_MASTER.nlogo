;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;TORTOISE RESTORATION TO
;;FLOREANA ISLAND
;;March 30, 2018
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;MODEL PURPOSE
;To simulate population dynamics of tortoise translocations while tracking population size in captivity and on islands, the percentage of native genomes captured, and total costs.

extensions [ r csv table ] ;need Windows to run r-netlogo extension, use csv extension to read in alleles, use table for creating look-up tables

;;;;;;;;;;;;;;;;;;;
;VARIABLES

globals[frac-female clutch-size-min clutch-size-max hatch-surv cost cost-release cost-Wolf corral adult-phi juv-phi young-phi release-eff-logit time-to-stop allele-data allele-data-loq allele-data-hiq counter num-microsat top5 mean-qval sd-qval aar sar ho n-fem n-male
           keepers-f adults.str adults.q ]

turtles-own[ sex age qval fecundity released? mom-id dad-id allele-info dad-alleles keep? maturity-age qval-eur ]

;;;;;;;;;;;;;;;;;;;;
;;SETUP PROCEDURES;;
;;;;;;;;;;;;;;;;;;;;

to setup
  ca

  set cost 35000 ;Cost for one "pre-adaptation" corral for the juveniles (max n=300) before they go to their islands. Start with one corral. If we reach more than 300 juveniles in captivity, need to build another.
  set corral 1
  set cost-release 10000 ;This is approx average between cost using barco Sierra Negra or La Molme
  set cost-Wolf 100000 ;Approximate cost of retrieving tortoises from Volcan Wolf (requiring helicopter support)

  set clutch-size-min 4
  set clutch-size-max 10

  set time-to-stop 0 ;This is used for the stopping criterion

  set num-microsat 21 ;number of microsatellite loci
  set mean-qval 0 ;Mean Q values only calculated at end of run
  set sd-qval 0

  ;Set survival rates: adult and juvenile survival rates drawn from Espanola distributions.  "Young" survival rate (in the wild) is between 0.6-0.9 (as used in Espanola paper)
  r:eval "load('C:/Users/ehunter/Dropbox/Galapagos/TortoiseRepatriation/Results/Espanola_PVA_2_short.RData') "
  set adult-phi r:get "Mod$sims.list$phiadult"
  set juv-phi r:get "Mod$sims.list$phijuv"
  set young-phi r:get "ys_samp"
  set frac-female r:get "ff_samp"
  set release-eff-logit r:get "Mod$sims.list$relbeta_phi"

  ;Get fertility distribution (product of egg viability and hatchling survival rate) from distribution created using approximate bayesian computation
  r:eval "load('C:/Users/ehunter/Dropbox/Galapagos/TortoiseRepatriation/Results/sims_ABC_fertoutput.RData') "
  set hatch-surv r:get "fertoutput"

  ;Set up "islands", create breeding corrals in captivity
  ask patches[
    if pxcor > 0 [set pcolor scale-color green pycor -2 20]  ;Pinta/Floreana (positive x coordinates)
    if pxcor < 0 [set pcolor scale-color blue pycor -2 20]   ;Santa Cruz breeding center (negative x)
  ]
  ask patches with [pycor = 20 and pxcor < -12]
    [set pcolor red]
  ask patches with [pycor = 20 and pxcor >= -12 and pxcor < -8]
    [set pcolor yellow]
  ask patches with [pycor = 20 and pxcor >= -8 and pxcor < -4]
    [set pcolor orange]
  ask patches with [pycor = 20 and pxcor >= -4 and pxcor < 0]
    [set pcolor brown]


  ;Read in allele values for Floreana hybrids, find number of breeders to create
  set allele-data csv:from-file "20_in_4Corrals_21_Loci_Genotypes.csv"
  set counter 0
  (foreach allele-data [
   crt 1[
   let temp-data item counter allele-data
   set counter counter + 1
   set sex item 1 temp-data ;"F" or "M"
   set qval item 2 temp-data
   set allele-info but-first but-first but-first temp-data
   ifelse sex = "F"
    [set color white
     set fecundity 0
     set mom-id first temp-data]
    [set color black
     set dad-id first temp-data]
   set age 20 ;unknown ages of breeders
   set maturity-age 20
   move-to one-of patches with [pcolor = blue]
   set ycor 20
   set keep? 1
   set dad-alleles [ ] ;make an empty list
   ]
  ])

  ;Have breeders move into their respective corrals
  ;Group 1=red, Group 2=yellow, Group 3=orange, Group 4=brown
  ;Females can go anywhere in their corrals, males are spaced out to help with corral-switching (their positioning is fixed right now, might want to make this random for skewed paternity scenarios)
  ;Group 1
  ask turtles with [mom-id = "Z7" or mom-id = "Z19" or mom-id = "Z9"] [move-to one-of patches with [pcolor = red]]
  ask turtles with [dad-id = "BR877" or dad-id = "I64"] [move-to one-of patches with [pcolor = red] set xcor -13]
  ;Group 2
  ask turtles with [mom-id = "BR106" or mom-id = "Z14" or mom-id = "Z21"] [move-to one-of patches with [pcolor = yellow]]
  ask turtles with [dad-id = "H210" or dad-id = "BR610"] [move-to one-of patches with [pcolor = yellow] set xcor -9]
  ;Group 3
  ask turtles with [mom-id = "I8" or mom-id = "Z11" or mom-id = "Z17"] [move-to one-of patches with [pcolor = orange]]
  ask turtles with [dad-id = "Z18" or dad-id = "I94"] [move-to one-of patches with [pcolor = orange] set xcor -5]
  ;Group 4
  ask turtles with [mom-id = "BR542" or mom-id = "E11" or mom-id = "G8"] [move-to one-of patches with [pcolor = brown]]
  ask turtles with [dad-id = "H2" or dad-id = "Z4"] [move-to one-of patches with [pcolor = brown] set xcor -1]

  ;Under scenario where high q-value adults are found on Wolf, increase corral capacity to 28.  Randomly choose 8 from 13 identified.  Randomly assign sexes (but assume 1 male, 1 female added to each corral).
  if(corral-capacity = 28)
  [ set cost cost + cost-Wolf ;Cost of retrieving tortoises from Wolf
    set allele-data-hiq csv:from-file "13_highQ_for_Additional_Breeders.csv"
    set allele-data-hiq n-of 8 allele-data-hiq  ;Randomly choose 8 of the 13
    set allele-data-hiq shuffle allele-data-hiq  ;Shuffle the list so that they go in different corrals each time
    set counter 0
    (foreach allele-data-hiq [
      crt 1[
      let temp-data item counter allele-data-hiq
      set counter counter + 1
      set allele-info but-first temp-data
      set age 20 ;unknown ages
      set maturity-age 20
      set ycor 20
      set keep? 1
      set dad-alleles [ ] ;make an empty list
      if(counter = 1)
      [ set sex "F"
        move-to one-of patches with [pcolor = red]]
      if(counter = 2)
      [ set sex "M"
        move-to one-of patches with [pcolor = red]]
      if(counter = 3)
      [ set sex "F"
        move-to one-of patches with [pcolor = yellow]]
      if(counter = 4)
      [ set sex "M"
        move-to one-of patches with [pcolor = yellow]]
      if(counter = 5)
      [ set sex "F"
        move-to one-of patches with [pcolor = orange]]
      if(counter = 6)
      [ set sex "M"
        move-to one-of patches with [pcolor = orange]]
      if(counter = 7)
      [ set sex "F"
        move-to one-of patches with [pcolor = brown]]
      if(counter = 8)
      [ set sex "M"
        move-to one-of patches with [pcolor = brown]]
      ifelse sex = "F"
      [set color white
       set fecundity 0
       set mom-id first temp-data] ;mom and dad IDs are not meaningful here, just for offspring
      [set color black
       set dad-id first temp-data] ]])
   ]

  reset-ticks
end



;;;;;;;;;;;;;;;;;
;;GO PROCEDURES;;
;;;;;;;;;;;;;;;;;


to go
  ;Breeding center reproduction only continues until release end, at end of releases, breeders are moved to Floreana
  if (ticks < release-end - release-age)
  [ask turtles with [ sex = "F" and age >= maturity-age and xcor < 0] [ reproduce-breeding-center ]]
  if (ticks >= release-end - release-age)
  [ask turtles with [ycor >= 20 and xcor < 0] [
      move-to one-of patches with [pcolor = green]
      set ycor 20]]

  ask turtles with [ sex = "F" and age >= maturity-age and xcor > 0] [ reproduce-wild ]

  ask turtles [ set released? 0 ] ;reset released variable for this year

  ask turtles [ survive ]

  ask turtles [
    set age age + 1
    if (age >= 1 and pycor < 20 ) [set ycor ycor + 1]]

  ;Move tortoises from captivity to Floreana
  if (any? turtles with [xcor < 0])
    [set cost cost + cost-release
      let releases turtles with [age > release-age and keep? = 0 and xcor < 0]
      set cost cost + (count releases * 15)  ;Genetic sampling occurs just before individuals are released
      ask releases [move-to one-of patches with [pcolor = green]
        set ycor age
        set released? 1]
    ]


  ;Releasing low q-value adults directly to Floreana
  if (release-adult-hybrids and ticks = 1) ;release lo q hybrids in year 1
  ;Read in allele values for Floreana hybrids, find number of breeders to create
  [set allele-data-loq csv:from-file "20_lowQ_for_Direct_Release.csv"
  set counter 0
  (foreach allele-data-loq [
   crt 1[
   let temp-data item counter allele-data-loq
   set counter counter + 1
   ;Sex unknown for these tortoises, so randomly assign each simulation
   set sex ifelse-value (random-float 1 >= 0.5) ["F"] ["M"]
   set allele-info but-first  temp-data
   ifelse sex = "F"
    [set color white
     set fecundity 0
     set mom-id first temp-data] ;mom and dad IDs are not meaningful here, just for offspring
    [set color black
     set dad-id first temp-data]
   set age 20 ;unknown ages at this point
   set maturity-age 20
   move-to one-of patches with [pcolor = green]
   set ycor 20
   set keep? 0
   set dad-alleles [ ] ;make an empty list
   ]
  ])
    set cost cost + cost-release + cost-Wolf
  ]

  ;Add in costs of captivity; if there are more than 300 juvenile tortoises in captivity, build a new pre-adaptation corral
  let captive-count count turtles with [xcor < 0 and age < 20]
  if (captive-count * corral) > (300 * corral)
  [set corral corral + 1
    set cost cost + 35000] ;Cost of new corral
  set cost cost + cost-captivity ;Cost of yearly tortoise care

  ;Add in a stop criterion
  if (stop-criterion-type = "time" and ticks >= 45) ;Time stopping criterion = 50 years
  [ set time-to-stop time-to-stop + 1 ]

  ;At end of run, output the alleles of all tortoises
  if (time-to-stop >= 6) ;Run the pop and cost options for 5 more years to let population dynamics settle out after final releases
  [ final-qvals
    stop]

  tick
end


;;;;;;;;;;;;;
;;SUBMODELS;;
;;;;;;;;;;;;;
to reproduce-wild
  ;select a mate
  let mate one-of turtles with [sex = "M" and age >= maturity-age and xcor > 0]
  if (mate = nobody) [ stop ]
  ;let genepool list ([genotype] of self) ([genotype] of mate)
  set dad-id [dad-id] of mate
  set dad-alleles [allele-info] of mate

  ;produce offspring, give genotype, and put in age class 1, use a poisson distribution to determine whether or not any recruitment
  let clutch-size random (clutch-size-max - clutch-size-min) + clutch-size-min
  let hatch-surv-draw one-of hatch-surv
  set fecundity random-poisson (hatch-surv-draw * clutch-size)

  if (fecundity >= 1) [
  let frac-female-draw one-of frac-female
  let n-females round(frac-female-draw * fecundity)
  let n-males fecundity - n-females
  if (n-females != nobody)
  [hatch n-females [
    move-to one-of patches with [pcolor = green]
    set ycor 0
    set age 0
    set maturity-age floor random-normal 25 2
    set fecundity 0
    set keep? 0
    set qval 0
    set qval-eur 0
    let mom-alleles allele-info
    set allele-info [ ]
    let n 0
    while [n < (num-microsat * 2)] [
      set allele-info lput (item (n + random 2) mom-alleles) allele-info ;Randomly select one of mom's alleles for the locus
      set allele-info lput (item (n + random 2) dad-alleles) allele-info ;dad's alleles
      set n n + 2
    ]
  ]]
  if (n-males != nobody)
  [hatch n-males [
    move-to one-of patches with [pcolor = green]
    set ycor 0
    set sex "M"
    set color black
    set age 0
    set maturity-age floor random-normal 25 2
    set keep? 0
    set qval 0
    set qval-eur 0
    let mom-alleles allele-info
    set allele-info [ ]
    let n 0
    while [n < (num-microsat * 2)] [
      set allele-info lput (item (n + random 2) mom-alleles) allele-info ;Randomly select one of mom's alleles for the locus
      set allele-info lput (item (n + random 2) dad-alleles) allele-info ;dad's alleles
      set n n + 2
    ]
  ]]]

  ;Clear dad-id at end
  set dad-id 0
end

to reproduce-breeding-center
  ;select a mate within my breeding corral (based on pcolor) and average genotypes
  let mate one-of turtles with [sex = "M" and pcolor = [pcolor] of myself]
  if (mate = nobody) [ stop ]
  set dad-id [dad-id] of mate
  set dad-alleles [allele-info] of mate

  ;produce offspring, give genotype, and put in age class 1 (need to determine rates still)
  let clutch-size random (clutch-size-max - clutch-size-min) + clutch-size-min
  let hatchlings clutch-size * 0.6 ;Eclosion rates of Espanola clutches in captivity is 0.6
  set fecundity hatchlings
  ;Sex-ratio can be controlled in captivity
  let n-females round(sex-ratio * fecundity)
  let n-males fecundity - n-females

  if (n-females != nobody) [
  hatch n-females [
    move-to one-of patches with [pcolor = blue]
    set ycor 0
    set age 0
    set maturity-age floor random-normal 25 2
    set keep? 0
    set qval 0
    set qval-eur 0
    let mom-alleles allele-info
    set allele-info [ ]
    let n 0
    while [n < (num-microsat * 2)] [
      set allele-info lput (item (n + random 2) mom-alleles) allele-info ;Randomly select one of mom's alleles for the locus
      set allele-info lput (item (n + random 2) dad-alleles) allele-info ;dad's alleles
      set n n + 2
    ]
  ]]
  if (n-males != nobody) [
  hatch n-males [
    move-to one-of patches with [pcolor = blue]
    set ycor 0
    set sex "M"
    set color black
    set age 0
    set maturity-age floor random-normal 25 2
    set keep? 0
    set qval 0
    set qval-eur 0
    let mom-alleles allele-info
    set allele-info [ ]
    let n 0
    while [n < (num-microsat * 2)] [
      set allele-info lput (item (n + random 2) mom-alleles) allele-info ;Randomly select one of mom's alleles for the locus
      set allele-info lput (item (n + random 2) dad-alleles) allele-info ;dad's alleles
      set n n + 2
    ]
  ]]

  ;Clear dad-id at end
  set dad-id 0
end


to survive
  if age < 5 [
    ifelse (xcor < 0)
    [let surv-prob 0.98 ;Survival of Espanola 1-5 year olds in captivity is ~98%
      if random-float 1 > surv-prob [die]]
    [let surv-prob-temp one-of young-phi
      let rel-eff-logit one-of release-eff-logit
      let rel-eff 1 / (1 + exp(-1 * rel-eff-logit))
      let surv-prob surv-prob-temp - (rel-eff * released?)
      if random-float 1 > surv-prob [die]]
    ]
  if age >= 5 and age < 8 [
    ifelse (xcor < 0)
    [let surv-prob 0.98
      if random-float 1 > surv-prob [die]]
    [let surv-prob-temp one-of juv-phi
      let rel-eff-logit one-of release-eff-logit
      let rel-eff 1 / (1 + exp(-1 * rel-eff-logit))
      let surv-prob surv-prob-temp - (rel-eff * released?)
    if random-float 1 > surv-prob [die]]
    ]
  if age >= 8 [
    ifelse (xcor < 0)
    [let surv-prob 0.995 ;adults in breeding center have very high survival
      if random-float 1 > surv-prob [die]]
    [let surv-prob one-of adult-phi
      if random-float 1 > surv-prob [die]]
    ]
end


to write-alleles
      let output (list self)
      foreach output
      [ ask ?
        [ file-write who
          file-type ","
          file-write mom-id
          file-type ","
          file-write dad-id
          file-type ","
          file-write keep?
          file-type ","
          file-write allele-info
          file-print ","]
      ]
end


;Calculates Q-values and diversity
to final-qvals
  let n-torts count turtles

  r:eval "library(ade4, lib.loc='C:/Users/ehunter/Documents/R/win-library/3.3')"
  r:eval "library(adegenet, lib.loc='C:/Users/ehunter/Documents/R/win-library/3.3')"
  r:eval "library(ape, lib.loc='C:/Users/ehunter/Documents/R/win-library/3.3')"
  r:eval "library(phangorn, lib.loc='C:/Users/ehunter/Documents/R/win-library/3.3')"
  r:eval "library(apex, lib.loc='C:/Users/ehunter/Documents/R/win-library/3.3')"
  r:eval "library(strataG, lib.loc='C:/Users/ehunter/Documents/R/win-library/3.3')"
  r:eval "setwd('C:/Users/ehunter/Dropbox/Galapagos/TortoiseRepatriation/Wolf/Data')"

  ;Read in the reference tortoises
  r:eval "ref <- read.csv('Reference_Data_STRUCTURE_REALLY_Reduced_ReBinned_NetLogo.csv')"
  r:eval "ref[ref==-9] <- NA"

  ;Append tortoise info to the reference tortoises (first delete old file)
  let file-name (word "C:/Users/ehunter/Desktop/TempFiles/temptorts" behaviorspace-run-number ".csv")
  file-open file-name
  ask turtles [write-alleles]
  file-close
  r:put "r_file_name" file-name
  r:eval "torts <- read.csv(r_file_name, header=FALSE)"
  r:eval "torts[torts==-9] <- NA"
  r:eval "alleles <- as.character(torts$V5)"
  r:eval "temp <- as.numeric(unlist(regmatches(alleles, gregexpr('[0-9]+|NA', alleles))))"
  r:eval "temp2 <- matrix(temp, nrow=nrow(torts), ncol=42, byrow=TRUE)"
  r:eval "temp2 <- as.data.frame(temp2)"
  r:eval "names(temp2) <- c('GAL45_1', 'GAL45_2', 'GAL50_1', 'GAL50_2', 'GAL75_1', 'GAL75_2', 'GAL94_1', 'GAL94_2', 'GAL100_1', 'GAL100_2', 'GAL127_1', 'GAL127_2', 'GAL136_1', 'GAL136_2', 'GAL159_1', 'GAL159_2', 'GAL263_1', 'GAL263_2', 'GAL194_1', 'GAL194_2', 'GAL288_1', 'GAL288_2', 'ACO63_1', 'ACO63_2', 'GAA45_1', 'GAA45_2', 'GAL158_1', 'GAL158_2', 'AC247_1', 'AC247_2', 'ACO39_1', 'ACO39_2', 'AC149_1', 'AC149_2', 'AC190_1', 'AC190_2', 'AGG68_1', 'AGG68_2', 'AC251_1', 'AC251_2', 'GAL21_1', 'GAL21_2')"
  r:eval "loci <- unique(substr(names(temp2), 1, 6))"

  ;#Append column name to allele numbers (except for NAs)
  r:eval "for (i in 1:ncol(temp2)){; for (j in 1:nrow(temp2)){; if(!is.na(temp2[j,i])) {temp2[j,i] <- paste(substr(names(temp2),1,6)[i], temp2[j,i], sep='_')}}}"

  ;Allelic richness:
  r:eval "aar <- length(unique(as.character(unlist(temp2))))"
  set aar r:get "aar"

  ;Shannon index of diversity (SAR)
  r:eval "alleles <- unique(as.vector(na.omit(unlist(temp2))))"
  r:eval "t3 <- as.vector(na.omit(unlist(temp2)))"
  r:eval "ps <- NULL"
  r:eval "ps2 <- NULL"
  r:eval "for(i in 1:length(alleles)){ps[i] <- length(t3[t3==alleles[i]]) / length(t3); ps2[i] <- ps[i] * log(ps[i])}"
  r:eval "sar <- -1 * sum(ps2)"
  set sar r:get "sar"

  ;To get other measures, needs to be in genpop format
  r:eval "library(ade4, lib.loc='C:/Users/ehunter/Documents/R/win-library/3.3')"
  r:eval "library(adegenet, lib.loc='C:/Users/ehunter/Documents/R/win-library/3.3')"
  r:eval "library(irelr, lib.loc='C:/Users/ehunter/Documents/R/win-library/3.3')"
  r:eval "library(data.table, lib.loc='C:/Users/ehunter/Documents/R/win-library/3.3')"
  r:eval "temp2 <- matrix(temp, nrow=nrow(torts), ncol=42, byrow=TRUE)"
  r:eval "cols <- seq(from=1, to=41, by=2)"
  r:eval "temp3 <- matrix(ncol=21, nrow=nrow(temp2))"
  r:eval "for (i in 1:(ncol(temp2)/2)) {temp3[,i] <- paste(temp2[,cols[i]], temp2[,cols[i]+1], sep='-')}"
  r:eval "flo2 <- df2genind(temp3, ploidy=2, sep='-')"

  ;Observed heterozygosity (Ho)
  r:eval "ho <- mean(summary(flo2)$Hobs)"
  set ho r:get "ho"

  ;Just take the first 12 loci to match with reference for STRUCTURE
  r:eval "temp2 <- temp2[,1:24]"
  r:eval "temp3 <- as.data.frame(cbind(torts[,1], rep(13, times=nrow(torts)), rep(0, times=nrow(torts)), temp2))"
  r:eval "temp3$V1 <- as.factor(temp3$V1)"
  r:eval "names(temp3) <- names(ref)"
  r:eval "ref_torts <- rbind(ref, temp3)"

  ;Run STRUCTURE
  r:eval "str.schemes <- as.data.frame(ref_torts[,2:3])"
  r:eval "rownames(str.schemes) <- ref_torts$Individual.ID"
  r:eval "msats <- new('gtypes', gen.data = ref_torts[,-(1:3)], ploidy=2, ind.names=as.vector(ref_torts[,1]), schemes=str.schemes, strata='Reference.Cluster')"
  r:eval "source('structureWriteNew2.R')"
  r:eval "sw.out <- structureWriteNew(msats, burnin=10000, numreps=50000, noadmix=FALSE, freqscorr=TRUE, label='TESTRUN/TESTRUN', popflag=ref_torts$POPFLAG)"
  r:eval "files <- sw.out$files"
  r:eval "exec <- 'structure'"
  r:eval "cmd <- paste(exec, ' -m ', files['mainparams'], ' -e ', files['extraparams'], ' -i ', files['data'], ' -o ', files['out'], sep = '')"
  r:eval "system(cmd)"

  ;Read in output and determine mean Q value
  r:eval "q.out <- read.table('C:/Users/ehunter/Dropbox/Galapagos/TortoiseRepatriation/Wolf/Data/TESTRUN/TESTRUN_out_q')"
  r:eval "q.out.fl <- q.out[99:109,3:9]"
  r:eval "idx <- which(q.out.fl>0.8, arr.ind=T)[1,2]" ;Determine which column corresponds to the Floreana cluster
  r:eval "q.out.hyb <- q.out[156:nrow(q.out),]"
  r:eval "hyb.q.fl <- q.out.hyb[,c(1,(idx+2))]"
  r:eval "names(hyb.q.fl) <- c('ID', 'Qval')"
  r:eval "qval.mn <- mean(hyb.q.fl$Qval)"
  r:eval "qval.sd <- sd(hyb.q.fl$Qval)"
  set mean-qval r:get "qval.mn"
  set sd-qval r:get "qval.sd"

end



;;;;;;;;;;;;;
;;REPORTERS;;
;;;;;;;;;;;;;


to-report pop-size-wild
  let wild-torts turtles with [xcor > 0]
  ifelse (count wild-torts = 0)
  [report 0]
  [let pop-size count wild-torts
  report pop-size]
end

to-report adult-pop-size-wild
  let wild-adults turtles with [xcor > 0 and age >= 15] ;For this purpose, larger tortoises have ecological impacts, so can use looser definition of "adult" (age=15)
  ifelse (count wild-adults = 0)
  [report 0]
  [let pop-size count wild-adults
  report pop-size]
end

to-report pop-size-captivity
  let captive-torts turtles with [xcor < 0]
  report count captive-torts
end

to-report cost-captivity ;For now, all tortoises (regardless of age) cost the same in captivity - average of $190/yr/tort
  let captive count turtles with [xcor < 0]
  let costs captive * 190 ; need to figure out cost per year per adult/juv
  report costs
end

to-report final-mean-qval
  report mean-qval
end

to-report final-sd-qval
  report sd-qval
end

to-report final-aar
  report aar
end

to-report final-sar
  report sar
end

to-report final-ho
  report ho
end

to-report final-sex-ratio
  let wild-torts turtles with [xcor > 0]
  let wild-fems wild-torts with [sex = "F"]
  report count wild-fems / count wild-torts
end






;
@#$#@#$#@
GRAPHICS-WINDOW
511
32
1181
483
16
-1
20.0
1
10
1
1
1
0
0
0
1
-16
16
0
20
0
0
1
ticks
30.0

BUTTON
123
97
187
130
Setup
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
223
97
286
130
Go
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
303
32
503
182
Population Size Pinta/Floreana
Time
Population Size
0.0
10.0
0.0
100.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot pop-size-wild"

PLOT
303
192
503
342
Costs
Time
Costs
0.0
10.0
0.0
500.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot cost"

SLIDER
120
153
292
186
release-end
release-end
10
100
50
10
1
NIL
HORIZONTAL

SWITCH
125
255
290
288
release-adult-hybrids
release-adult-hybrids
0
1
-1000

PLOT
303
355
503
505
Population Size Captivity
Time
Population Size
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot pop-size-captivity"

CHOOSER
154
332
292
377
stop-criterion-type
stop-criterion-type
"time" "pop" "cost"
0

SLIDER
120
187
292
220
release-age
release-age
3
15
3
1
1
NIL
HORIZONTAL

SLIDER
120
221
292
254
sex-ratio
sex-ratio
0
1
0.5
0.05
1
NIL
HORIZONTAL

CHOOSER
154
379
292
424
corral-capacity
corral-capacity
20 28
1

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270

@#$#@#$#@
NetLogo 5.3.1
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="AllOptions" repetitions="5" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>pop-size-wild</metric>
    <metric>adult-pop-size-wild</metric>
    <metric>cost</metric>
    <metric>final-mean-qval</metric>
    <metric>final-sd-qval</metric>
    <metric>final-aar</metric>
    <metric>final-sar</metric>
    <metric>final-ho</metric>
    <metric>final-sex-ratio</metric>
    <metric>ticks</metric>
    <enumeratedValueSet variable="release-adult-hybrids">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-criterion-type">
      <value value="&quot;time&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="release-end">
      <value value="20"/>
      <value value="30"/>
      <value value="40"/>
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="release-age">
      <value value="3"/>
      <value value="4"/>
      <value value="5"/>
      <value value="6"/>
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sex-ratio">
      <value value="0.5"/>
      <value value="0.67"/>
      <value value="0.75"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Corral-capacity">
      <value value="20"/>
      <value value="28"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Optimals" repetitions="20" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>pop-size-wild</metric>
    <metric>adult-pop-size-wild</metric>
    <metric>cost</metric>
    <metric>final-aar</metric>
    <metric>final-sar</metric>
    <metric>final-ho</metric>
    <metric>final-sex-ratio</metric>
    <metric>ticks</metric>
    <enumeratedValueSet variable="release-adult-hybrids">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-criterion-type">
      <value value="&quot;time&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="release-end">
      <value value="20"/>
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="release-age">
      <value value="3"/>
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sex-ratio">
      <value value="0.67"/>
      <value value="0.75"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Corral-capacity">
      <value value="28"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="20_rel_end" repetitions="5" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>pop-size-wild</metric>
    <metric>adult-pop-size-wild</metric>
    <metric>cost</metric>
    <metric>final-aar</metric>
    <metric>final-sar</metric>
    <metric>final-ho</metric>
    <metric>final-sex-ratio</metric>
    <metric>ticks</metric>
    <enumeratedValueSet variable="release-adult-hybrids">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-criterion-type">
      <value value="&quot;time&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="release-end">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="release-age">
      <value value="3"/>
      <value value="4"/>
      <value value="5"/>
      <value value="6"/>
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sex-ratio">
      <value value="0.5"/>
      <value value="0.67"/>
      <value value="0.75"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Corral-capacity">
      <value value="20"/>
      <value value="28"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="AllOptions2" repetitions="5" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>pop-size-wild</metric>
    <metric>adult-pop-size-wild</metric>
    <metric>cost</metric>
    <metric>final-mean-qval</metric>
    <metric>final-sd-qval</metric>
    <metric>final-aar</metric>
    <metric>final-sar</metric>
    <metric>final-ho</metric>
    <metric>final-sex-ratio</metric>
    <metric>ticks</metric>
    <enumeratedValueSet variable="release-adult-hybrids">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-criterion-type">
      <value value="&quot;time&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="release-end">
      <value value="40"/>
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="release-age">
      <value value="3"/>
      <value value="4"/>
      <value value="5"/>
      <value value="6"/>
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sex-ratio">
      <value value="0.5"/>
      <value value="0.67"/>
      <value value="0.75"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Corral-capacity">
      <value value="20"/>
      <value value="28"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Optimals_diversity" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>pop-size-wild</metric>
    <metric>adult-pop-size-wild</metric>
    <metric>cost</metric>
    <metric>final-mean-qval</metric>
    <metric>final-sd-qval</metric>
    <metric>final-aar</metric>
    <metric>final-sar</metric>
    <metric>final-ho</metric>
    <metric>final-sex-ratio</metric>
    <metric>ticks</metric>
    <enumeratedValueSet variable="release-adult-hybrids">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stop-criterion-type">
      <value value="&quot;time&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="release-end">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="release-age">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sex-ratio">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Corral-capacity">
      <value value="28"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180

@#$#@#$#@
0
@#$#@#$#@
