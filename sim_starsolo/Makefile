include Mf.common

.SECONDARY:

.PHONY: all

all: real sims_run


real:
	$(MAKE) All_Gene_Real       genome=human_CR_3.0.0   sample=10X/3/pbmc_5k   threadRun=20/run1
	# TODO:	$(MAKE) All_GeneFull_Real   genome=human_CR_3.0.0   sample=10X/3/pbmc_5k   threadRun=20/run1

sims_run:
	$(MAKE) All_Gene_Sims       genome=human_CR_3.0.0   sample=10X/3/pbmc_5k_sims_human_CR_3.0.0_MultiGeneNo    threadRun=20/run1
	# TODO: $(MAKE) All_Gene_Sims       genome=human_CR_3.0.0   sample=10X/3/pbmc_5k_sims_human_CR_3.0.0_MultiGeneYes   threadRun=20/run1
	# TODO: $(MAKE) All_Gene_Sims       genome=human_CR_3.0.0   sample=10X/3/pbmc_5k_sims_human_CR_3.0.0_MultiGeneNo_OnlyExonicReads   threadRun=20/run1

############################################################ Versions
alevin-fry  := salmon-alevin-fry_1.9.0_0.8.0
kbpy        := kallisto_0.49.0
STAR        := STAR_2.7.9a
CellRanger3 := CellRanger_3.1.0
CellRanger6 := CellRanger_6.1.2
CellRanger7 := CellRanger_7.0.1

############################################################ tools and options
# Gene: standard gene expression

All_GeneFull_Real:   \
                     GeneFull/STAR_CR4 \
                     GeneFull/CellRanger6 \
                     GeneFull/CellRanger7
	
All_Gene_Sims: Gene/kbpy \
               Gene/STAR_sims \
               #TODO: Gene/alevin-fry

All_Gene_Real:  \
                Gene/kbpy \
		Gene/alevin-fry \
		#TODO :Gene/STAR_CR4 \
		#TODO: Gene/CellRanger6 \
		#TODO: Gene/CellRanger7

Gene/CellRanger6: 
	make -f Mf.count count/$(CellRanger6)/$(genome)/standard/default/$(sample)/$(threadRun)/CO

GeneFull/CellRanger6: 
	make -f Mf.count count/$(CellRanger6)/$(genome)/standard/intronic/$(sample)/$(threadRun)/CO

Gene/CellRanger7: 
	make -f Mf.count count/$(CellRanger7)/$(genome)/standard/default/$(sample)/$(threadRun)/CO

GeneFull/CellRanger7: 
	make -f Mf.count count/$(CellRanger7)/$(genome)/standard/intronic/$(sample)/$(threadRun)/CO

Gene/STAR_sims:
	$(MAKE) -f Mf.count count/$(STAR)/$(genome)/fullSA/10X_noSAM_sims_mult_ENCODE/$(sample)/$(threadRun)/CO

Gene/STAR_CR4:
	make -f Mf.count count/$(STAR)/$(genome)/fullSA/10X_CR4_noSAM/$(sample)/$(threadRun)/CO
	make -f Mf.count count/$(STAR)/$(genome)/sparseSA3/10X_CR4_noSAM/$(sample)/$(threadRun)/CO

GeneFull/STAR_CR4:
	$(MAKE) -f Mf.count count/$(STAR)/$(genome)/fullSA/10X_CR4_GeneFull_noSAM/$(sample)/$(threadRun)/CO
	$(MAKE) -f Mf.count count/$(STAR)/$(genome)/sparseSA3/10X_CR4_GeneFull_noSAM/$(sample)/$(threadRun)/CO

Gene/kbpy:
	make -f Mf.count count/$(kbpy)/$(genome)/standard_1/default/$(sample)/$(threadRun)/CO
	make -f Mf.count count/$(kbpy)/$(genome)/standard_offlist_1/default/$(sample)/$(threadRun)/CO
	#TODO make -f Mf.count count/$(kbpy)/$(genome)/standard_1/mult/$(sample)/$(threadRun)/CO
	#TODO :make -f Mf.count count/$(kbpy)/$(genome)/nucleus_1/default/$(sample)/$(threadRun)/CO
	# TODO: Add more kallisto configurations here

Gene/alevin-fry:
	make -f Mf.count count/$(alevin-fry)/$(genome)/standard/sketch_rad/$(sample)/$(threadRun)/CO afGplOpt=knee
	make -f Mf.count count/$(alevin-fry)/$(genome)/standard/rad/$(sample)/$(threadRun)/CO        afGplOpt=knee
	# make -f Mf.count count/$(alevin-fry)/$(genome)/decoyFull/rad/$(sample)/$(threadRun)/CO       afGplOpt=knee
	# make -f Mf.count count/$(alevin-fry)/$(genome)/decoyPartial/rad/$(sample)/$(threadRun)/CO    afGplOpt=knee
	make -f Mf.count count/$(alevin-fry)/$(genome)/splici/rad/$(sample)/$(threadRun)/CO          afGplOpt=knee
	make -f Mf.count count/$(alevin-fry)/$(genome)/splici/sketch_rad/$(sample)/$(threadRun)/CO          afGplOpt=knee
	make -f Mf.count count/$(alevin-fry)/$(genome)/splici_sparse/rad/$(sample)/$(threadRun)/CO          afGplOpt=knee
	make -f Mf.count count/$(alevin-fry)/$(genome)/splici_sparse/sketch_rad/$(sample)/$(threadRun)/CO          afGplOpt=knee


