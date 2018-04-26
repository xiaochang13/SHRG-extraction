#!/bin/bash
#python -u ./hrg_learn.py --output_dir rule_dir --data_dir ./simple_train --stats_dir stats --use_stats
python -u ./hrg_learn.py --output_dir rule_dir --data_dir ./train --stats_dir stats --use_stats --format amr > amr_labeled.log &
#python -u ./hrg_learn.py --output_dir rule_dir --data_dir ./train --stats_dir stats --use_stats --random
# python -u ./hrg_learn.py --output_dir rule_dir --data_dir ./train --stats_dir stats --use_stats --nolabel > direction_only.log & 
# python -u ./hrg_learn.py --output_dir rule_dir --data_dir ./train --stats_dir stats --use_stats --nolabel --nodir > unlabeled.log & 
# python -u ./hrg_learn.py --output_dir rule_dir --data_dir ./train --stats_dir stats --use_stats --format amr > labeled.log &

# Dependency results
# python -u ./hrg_learn.py --output_dir dep_output/english_test --conll_file ./dependency/conll_2007_ara_eng/data/english/ptb/train/english_ptb_train.conll --stats_dir stats --use_stats --format dep
# python -u ./hrg_learn.py --output_dir dep_output/english --conll_file ./dependency/conll_2007_ara_eng/data/english/ptb/train/english_ptb_train.conll --stats_dir stats --use_stats --format dep > english.log &
# python -u ./hrg_learn.py --output_dir dep_output/arabic --conll_file ./dependency/conll_2007_ara_eng/data/arabic/PADT/train/arabic_PADT_train.conll --stats_dir stats --use_stats --format dep > arabic.log &
python -u ./hrg_learn.py --output_dir dep_output/italian --conll_file ./dependency/conll_2007_gre_hung_ita/data/italian/isst/train/italian_isst_train.conll --stats_dir stats --use_stats --format dep > italian.log &
# python -u ./hrg_learn.py --output_dir dep_output/czech --conll_file ./dependency/conll_2007_bas_cat_tur_cze/data/czech/pdt/train/czech_pdt_train.conll --stats_dir stats --use_stats --format dep > czech.log &
# python -u ./hrg_learn.py --output_dir dep_output/czech --conll_file ./dependency/conll_2007_bas_cat_tur_cze/data/czech/pdt/train/czech_pdt_train.conll --stats_dir stats --use_stats --format dep
