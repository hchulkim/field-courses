#####################################################################
# Title: Makefile for analysis report
# Maintainer: Hyoungchul
# Initial date: 2025-11-13
# Modified date: 2025-11-13
#####################################################################

## Directory vars (usually only these need changing)
rdir = code/hyoungchul/
datadir = data/raw/
outputdir = output/
figdir    = output/figures/
tabdir    = output/tables/

.PHONY: all create

# make all execution formula
all: create

## Draw the Makefile DAG
## Requires: https://github.com/lindenb/makefile2graph
dag: makefile-dag.png

makefile-dag.png: Makefile
	@$(MAKE) -Bnd all | make2graph | dot -Tpng -Gdpi=300 -o $@

#####################################################################
# running R codes
#####################################################################

create: $(tabdir)author_table_c1_a.tex $(figdir)author_fig2.png $(tabdir)author_table3.tex $(tabdir)author_table4.tex $(tabdir)author_table5.tex $(figdir)author_appendix_fig_a2.png $(tabdir)author_appendix_table_b1.tex $(tabdir)author_appendix_table_c11.tex $(tabdir)author_appendix_table_c11_200.tex $(figdir)grid_count_histogram.png

$(tabdir)author_table_c1_a.tex $(figdir)author_fig2.png $(tabdir)author_table3.tex $(tabdir)author_table4.tex $(tabdir)author_table5.tex &: $(datadir)PanelSouthChina.dta $(datadir)SouthChinaPixelLevelMain.dta $(datadir)CountyLevelDataset.dta $(rdir)/rerun_original_main_code.R
	Rscript $(rdir)rerun_original_main_code.R

$(figdir)author_appendix_fig_a2.png $(tabdir)author_appendix_table_b1.tex $(tabdir)author_appendix_table_c11.tex $(tabdir)author_appendix_table_c11_200.tex &: $(datadir)MannEtAl2009Data.dta $(datadir)CountyLevelDataset.dta $(datadir)SouthChinaPixelLevelMain.dta $(rdir)rerun_original_appendix_code.R
	Rscript $(rdir)rerun_original_appendix_code.R

$(figdir)grid_count_histogram.png : $(datadir)PanelSouthChina.dta $(rdir)explore_seats_data.R
	Rscript $(rdir)explore_seats_data.R
