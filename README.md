## Группа G4 / SRR11461739, SRR11461738
Участники:  
* Соколов Василий Викторович
* Терентьева Юлия Андреевна
* Фаттахов Айдар Ильгизович
* Федоров Александр Николаевич

### Шаг 0 - подготовка окружения
Работать будем с пакетными менеджером **conda**, для чего-то другого(например Docker-a) все равно прав нет. Ставим все искомые утилиты:
```shell
# Новое окружение
conda create -y --name nanopore-hw
conda activate nanopore-hw

# Пакеты
conda install -y -c bioconda \
	sra-tools \ # fastq-dump 
	fastqc    \ # qc
	multiqc   \ # qc
	porechop  \ # nanopore trimming
	filtlong  \ # long reads filtering
	unicycler \ # bacterial genome assembler
	abricate    # screening for antimicrobial resistance or virulence genes
```

К сожалению, bioconda рецепт для r-minionqc(qc для длинных прочтений) сломан, r-minionqc не ставится даже в чистое окружение. 
Поставим зависимости через conda, а скрипт просто скачаем:

```shell
conda install -y r-data.table r-futile.logger r-ggplot2 r-optparse r-plyr r-readr r-reshape2 r-scales r-viridis r-yaml
wget https://raw.githubusercontent.com/roblanf/minion_qc/master/MinIONQC.R -O MinIONQC.R
```
### Шаг 1 - загрузка прочтений
Выполним обязательно настройку для `fastq-dump`. В нашем случае это просто save-exit, не очень понятно зачем она сделали ее обязательной.

```shell
vdb-config --interactive
```

Загружаем наши прочтения. Использовать gzip не будем, на кластере свободно порядка 3TB:
```shell
mkdir -p fastq/SRR11461738 fastq/SRR11461739

# paired end
fastq-dump --split-files -O fastq/SRR11461738 SRR11461738

# long reads
fastq-dump -O fastq/SRR11461739 SRR11461739
```
### Шаг 3 - предобработка

Сделаем подвыборку длиных прочтений (--length_weight 10 для приоритета именно длинных прочтений):
```shell
filtlong --min_length 1000 --target_bases 400000000 --length_weight 10 fastq/SRR11461739/SRR11461739.fastq > fastq/SRR11461739/SRR11461739.filtlong.fastq
```

Если судить по [twitter](https://twitter.com/rrwick/status/1134058857334886401?lang=en)(нормальной документации в google не нашлось) 
то guppy уже это делает автоматически. Наши прочтения с 2020-го года, можно считать, что адаптеров там нет. 
Однако, на всякий случай, воспользуемся porechop
```shell
porechop -i fastq/SRR11461739/SRR11461739.filtlong.fastq -o fastq/SRR11461739/SRR11461739.trimmed.filtlong.fastq
```

Забавно, адаптеры нашлись, притом почти везде:
```shell
27,031 / 29,786 reads had adapters trimmed from their start (1,001,422 bp removed)
24,766 / 29,786 reads had adapters trimmed from their end (371,069 bp removed)
3,291 / 29,786 reads were split based on middle adapters
```

В идеале porechop надо запускать на исходных прочтениях, и только потом делать filtlong. Однако для такого подхода 
на сервере не хватает RAM, а сделать swap нет прав. Оставим как есть, но в боевой ситуации нужно было бы что-то придумать
(купить RAM/написать сис. админу/разбить файл на несколько/т.д.).

### Шаг 4 - qc

Сначала fastqc для парных прочтений:
```shell
mkdir -p qc/fastqc qc/minionqc

fastqc -t 16 -o qc/fastqc fastq/SRR11461738/SRR11461738_* fastq/SRR11461739/SRR11461739.*
```

Попробуем minionqc на исходных/отфильтрованных/обрезанных прочтениях:
```shell
Rscript MinIONQC.R -p 16 -i fastq/SRR11461739/SRR11461739.fastq -o qc/minionqc/SRR11461739
Rscript MinIONQC.R -p 16 -i fastq/SRR11461739/SRR11461739.filtlong.fastq -o qc/minionqc/SRR11461739.filtlong
Rscript MinIONQC.R -p 16 -i fastq/SRR11461739/SRR11461739.trimmed.filtlong.fastq -o qc/minionqc/SRR11461739.trimmed.filtlong.fastq
```
К сожалению, во всех случаях minionqc упал с ошибкой следующего вида:
```shell
INFO [2020-11-03 22:17:20] Loading input file: fastq/SRR11461739/SRR11461739.filtlong.fastq
WARN [2020-11-03 22:17:26] There were problems in parsing your input file with the following rows: 
WARN [2020-11-03 22:17:26] 895
WARN [2020-11-03 22:17:26] 895
WARN [2020-11-03 22:17:26] 895
# Повтор одного и того же сообщения тысячи раз...
WARN [2020-11-03 22:17:26] 1095
# Повтор одного и того же сообщения тысячи раз... + еще десяток подобных ошибок
Error: Assigned data `as.numeric(as.character(d$sequence_length_template))` must be compatible with existing data.
✖ Existing data has 64525 rows.
✖ Assigned data has 0 rows
Backtrace:
     █
  1. └─global::single.flowcell(input.file, output.dir, q)
  2.   └─global::load_summary(input.file, min.q = c(-Inf, q))
  3.     ├─base::`$<-`(`*tmp*`, "sequence_length_template", value = numeric(0))
  4.     └─tibble:::`$<-.tbl_df`(`*tmp*`, "sequence_length_template", value = numeric(0))
  5.       └─tibble:::tbl_subassign(...)
  6.         └─tibble:::vectbl_recycle_rhs(...)
  7.           └─base::tryCatch(...)
  8.             └─base:::tryCatchList(expr, classes, parentenv, handlers)
  9.               └─base:::tryCatchOne(expr, names, parentenv, handlers[[1L]])
 10.        
Execution halted
```
В github репозитории висит [issue](https://github.com/roblanf/minion_qc/issues/49) на подобную ошибку с пометкой исправлено. 
Только наша ошибка - это падение в уже исправленном коде. Т.к. это учебный проект, поверим что прочтения хорошие и углубляться 
в проблему не будем.  
*P.S. Казалось бы, есть ли связь между языком и качеством утилит? А ведь тулы на R очень часто из коробки не работают...*

Агрегируем qc-отчет:
```shell
multiqc -o qc qc
```
MultiQC отчет доступен по [ссылке](https://htmlpreview.github.io/?https://github.com/alnfedorov/nanopore-hw/blob/master/qc/multiqc_report.html).

Из отчета можно заметить, что адаптеры у парных прочтений уже обрезаны, однако есть странное распределение в первых 20 и последних 5 позициях. 
Кажется, что это просто артефакт секвенирования на Illumina(будет здорово, если Вы прокомментируете). 
Чтобы не потерять информацию, обрезать их не будем. Если сборка получится плохой, то всегда можно вернуться и попробовать обрезать.

Для предобработанных длинных прочтений такое почти не характерно, но 5\` конец также можно обрезать при необходимости.  
(Странное распределение на 3\` конце можно объяснить очень просто - настолько длинных прочтений всего несколько.)
