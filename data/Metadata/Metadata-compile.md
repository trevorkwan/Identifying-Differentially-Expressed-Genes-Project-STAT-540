Metadat\_compile
================
Yilin Qiu
15/03/2021

#### download dataset from the Gene Expression Omnibus (GEO)

``` r
# Download info for each sample 
metadata <- getGEO(filename = 'GSE157344_family.soft' )
```

    ## Reading file....

    ## Parsing....

    ## Found 55 entities...

    ## GPL18573 (1 of 56 entities)

    ## GSM4762139 (2 of 56 entities)

    ## GSM4762140 (3 of 56 entities)

    ## GSM4762141 (4 of 56 entities)

    ## GSM4762142 (5 of 56 entities)

    ## GSM4762143 (6 of 56 entities)

    ## GSM4762144 (7 of 56 entities)

    ## GSM4762145 (8 of 56 entities)

    ## GSM4762146 (9 of 56 entities)

    ## GSM4762147 (10 of 56 entities)

    ## GSM4762148 (11 of 56 entities)

    ## GSM4762149 (12 of 56 entities)

    ## GSM4762150 (13 of 56 entities)

    ## GSM4762151 (14 of 56 entities)

    ## GSM4762152 (15 of 56 entities)

    ## GSM4762153 (16 of 56 entities)

    ## GSM4762155 (17 of 56 entities)

    ## GSM4762156 (18 of 56 entities)

    ## GSM4762157 (19 of 56 entities)

    ## GSM4762158 (20 of 56 entities)

    ## GSM4762159 (21 of 56 entities)

    ## GSM4762160 (22 of 56 entities)

    ## GSM4762161 (23 of 56 entities)

    ## GSM4762162 (24 of 56 entities)

    ## GSM4762163 (25 of 56 entities)

    ## GSM4762164 (26 of 56 entities)

    ## GSM4762165 (27 of 56 entities)

    ## GSM4762166 (28 of 56 entities)

    ## GSM4762167 (29 of 56 entities)

    ## GSM4762168 (30 of 56 entities)

    ## GSM4762169 (31 of 56 entities)

    ## GSM4762170 (32 of 56 entities)

    ## GSM4762171 (33 of 56 entities)

    ## GSM4762172 (34 of 56 entities)

    ## GSM4762173 (35 of 56 entities)

    ## GSM4762174 (36 of 56 entities)

    ## GSM4762175 (37 of 56 entities)

    ## GSM4762177 (38 of 56 entities)

    ## GSM4762178 (39 of 56 entities)

    ## GSM4762179 (40 of 56 entities)

    ## GSM4762180 (41 of 56 entities)

    ## GSM4762181 (42 of 56 entities)

    ## GSM4762182 (43 of 56 entities)

    ## GSM4762183 (44 of 56 entities)

    ## GSM4762184 (45 of 56 entities)

    ## GSM4762185 (46 of 56 entities)

    ## GSM4762186 (47 of 56 entities)

    ## GSM4762187 (48 of 56 entities)

    ## GSM4762188 (49 of 56 entities)

    ## GSM4762189 (50 of 56 entities)

    ## GSM4762190 (51 of 56 entities)

    ## GSM4762191 (52 of 56 entities)

    ## GSM4762192 (53 of 56 entities)

    ## GSM4762193 (54 of 56 entities)

    ## GSM4762194 (55 of 56 entities)

``` r
# Select the first sample and extract the characteristics for the first sample
metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762139"]]@header[["characteristics_ch1"]])
# Split one column  into two columns
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
# Convert columns into rows
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
# Add Sample_id 
metadata_3$sample = metadata@gsms[["GSM4762139"]]@header[["geo_accession"]]
# Join the first sample to a black dataframe
metadata_coml <- rbind(metadata_3)


#####Select the second sample and extract the characteristics for the second sample. This will be repeated until the last sample has been extracted. 
metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762140"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762140"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3, metadata_coml )



metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762141"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762141"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3, metadata_coml )


metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762142"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762142"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)



metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762143"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762143"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)



metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762144"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762144"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)




metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762145"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762145"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)


metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762146"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762146"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)


metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762147"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762147"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)


metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762148"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762148"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)
```

``` r
metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762149"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762149"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)


metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762150"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762150"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)


metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762151"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762151"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)


metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762152"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762152"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)


metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762153"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762153"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)



metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762155"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762155"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)


metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762156"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762156"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)


metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762157"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762157"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)


metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762158"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762158"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)



metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762159"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762159"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)
```

``` r
metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762160"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762160"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)


metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762161"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762161"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)



metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762162"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762162"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)



metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762163"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762163"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)



metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762164"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762164"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)



metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762165"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762165"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)



metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762166"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762166"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)



metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762167"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762167"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)


metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762168"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762168"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)



metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762169"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762169"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)
```

``` r
metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762170"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762170"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)


metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762171"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762171"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)



metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762172"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762172"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)



metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762173"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762173"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)



metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762174"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762174"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)



metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762175"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762175"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)





metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762177"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762177"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)


metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762178"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762178"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)



metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762179"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762179"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)


metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762180"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762180"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)
```

``` r
metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762181"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762181"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)


metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762182"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762182"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)



metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762183"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762183"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)



metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762184"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762184"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)



metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762185"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762185"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)



metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762186"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762186"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)





metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762187"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762187"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)


metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762188"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762188"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)



metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762189"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762189"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)


metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762190"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762190"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)




metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762191"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762191"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)




metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762192"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762192"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)




metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762193"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762193"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)




metadata_2 <- data.frame(characteristic = metadata@gsms[["GSM4762194"]]@header[["characteristics_ch1"]])
metadata_2 <- metadata_2 %>% separate(characteristic, c("x", "y"), ": ")
metadata_3<-as.data.frame(t(metadata_2)) 
names(metadata_3)<- metadata_3[1,]
metadata_3 <- metadata_3[-1,]
metadata_3$sample = metadata@gsms[["GSM4762194"]]@header[["geo_accession"]]
metadata_coml <- rbind(metadata_3,  metadata_coml)


metadata_comlete <- rbind(metadata_coml)
str(metadata_comlete)
```

    ## 'data.frame':    54 obs. of  8 variables:
    ##  $ clinic status   : chr  "Severe COVID" "Severe COVID" "Severe COVID" "Severe COVID" ...
    ##  $ tissue          : chr  "Blood" "Blood" "Blood" "Blood" ...
    ##  $ subject id      : chr  "Patient 9" "Patient 8" "Patient 7" "Patient 6" ...
    ##  $ sofa score      : chr  "5" "4" "5" "6" ...
    ##  $ age             : chr  "58" "55" "70" "59" ...
    ##  $ Sex             : chr  "Male" "Male" "Male" "Female" ...
    ##  $ clinical outcome: chr  "Alive" "Alive" "Dead" "Alive" ...
    ##  $ sample          : chr  "GSM4762194" "GSM4762193" "GSM4762192" "GSM4762191" ...

``` r
write.csv(metadata_comlete, "metadata_comlete.csv")
```
