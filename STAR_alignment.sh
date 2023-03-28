#!/bin/bash

STAR --runThreadN 12 --genomeDir /mnt/c/RH22023261/index_mm10 --readFilesIn /mnt/c/RH22023261/data/AMT_E_1_1.fq /mnt/c/RH22023261/data/AMT_E_1_2.fq --outFileNamePrefix AMT_E_1 --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts

STAR --runThreadN 12 --genomeDir /mnt/c/RH22023261/index_mm10 --readFilesIn /mnt/c/RH22023261/data/AMT_E_2_1.fq /mnt/c/RH22023261/data/AMT_E_2_2.fq --outFileNamePrefix AMT_E_2 --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts

STAR --runThreadN 12 --genomeDir /mnt/c/RH22023261/index_mm10 --readFilesIn /mnt/c/RH22023261/data/MT_E_1_1.fq /mnt/c/RH22023261/data/MT_E_1_2.fq --outFileNamePrefix MT_E_1 --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts

STAR --runThreadN 12 --genomeDir /mnt/c/RH22023261/index_mm10 --readFilesIn /mnt/c/RH22023261/data/MT_E_2_1.fq /mnt/c/RH22023261/data/MT_E_2_2.fq --outFileNamePrefix MT_E_2 --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts

STAR --runThreadN 12 --genomeDir /mnt/c/RH22023261/index_mm10 --readFilesIn /mnt/c/RH22023261/data/AMT_N_1_1.fq /mnt/c/RH22023261/data/AMT_N_1_2.fq --outFileNamePrefix AMT_N_1 --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts

STAR --runThreadN 12 --genomeDir /mnt/c/RH22023261/index_mm10 --readFilesIn /mnt/c/RH22023261/data/AMT_N_2_1.fq /mnt/c/RH22023261/data/AMT_N_2_2.fq --outFileNamePrefix AMT_N_2 --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts

STAR --runThreadN 12 --genomeDir /mnt/c/RH22023261/index_mm10 --readFilesIn /mnt/c/RH22023261/data/MT_N_1_1.fq /mnt/c/RH22023261/data/MT_N_1_2.fq --outFileNamePrefix MT_N_1 --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts

STAR --runThreadN 12 --genomeDir /mnt/c/RH22023261/index_mm10 --readFilesIn /mnt/c/RH22023261/data/MT_N_2_1.fq /mnt/c/RH22023261/data/MT_N_2_2.fq --outFileNamePrefix MT_N_2 --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts

STAR --runThreadN 12 --genomeDir /mnt/c/RH22023261/index_mm10 --readFilesIn /mnt/c/RH22023261/data/VIVO_E_1_1.fq /mnt/c/RH22023261/data/VIVO_E_1_2.fq --outFileNamePrefix VIVO_E_1 --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts

STAR --runThreadN 12 --genomeDir /mnt/c/RH22023261/index_mm10 --readFilesIn /mnt/c/RH22023261/data/VIVO_E_2_1.fq /mnt/c/RH22023261/data/VIVO_E_2_2.fq --outFileNamePrefix VIVO_E_2 --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts

STAR --runThreadN 12 --genomeDir /mnt/c/RH22023261/index_mm10 --readFilesIn /mnt/c/RH22023261/data/VIVO_E_3_1.fq /mnt/c/RH22023261/data/VIVO_E_3_2.fq --outFileNamePrefix VIVO_E_3 --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts

STAR --runThreadN 12 --genomeDir /mnt/c/RH22023261/index_mm10 --readFilesIn /mnt/c/RH22023261/data/VIVO_E_4_1.fq /mnt/c/RH22023261/data/VIVO_E_4_2.fq --outFileNamePrefix VIVO_E_4 --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts

STAR --runThreadN 12 --genomeDir /mnt/c/RH22023261/index_mm10 --readFilesIn /mnt/c/RH22023261/data/VIVO_N_1_1.fq /mnt/c/RH22023261/data/VIVO_N_1_2.fq --outFileNamePrefix VIVO_N_1 --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts

STAR --runThreadN 12 --genomeDir /mnt/c/RH22023261/index_mm10 --readFilesIn /mnt/c/RH22023261/data/VIVO_N_2_1.fq /mnt/c/RH22023261/data/VIVO_N_2_2.fq --outFileNamePrefix VIVO_N_2 --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts

STAR --runThreadN 12 --genomeDir /mnt/c/RH22023261/index_mm10 --readFilesIn /mnt/c/RH22023261/data/VIVO_N_3_1.fq /mnt/c/RH22023261/data/VIVO_N_3_2.fq --outFileNamePrefix VIVO_N_3 --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts

STAR --runThreadN 12 --genomeDir /mnt/c/RH22023261/index_mm10 --readFilesIn /mnt/c/RH22023261/data/VIVO_N_4_1.fq /mnt/c/RH22023261/data/VIVO_N_4_2.fq --outFileNamePrefix VIVO_N_4 --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts