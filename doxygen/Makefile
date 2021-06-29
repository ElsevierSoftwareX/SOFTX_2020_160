
.PHONY: filter folders images doxygen

export DOXYGEN_OUTPUT_DIRECTORY := ../public
export DOXYGEN_SOURCE_DIRECTORY := ./
export DOXYGEN_EXAMPLE_PATH := $(DOXYGEN_OUTPUT_DIRECTORY)/macros
export DOXYGEN_IMAGE_PATH := $(DOXYGEN_OUTPUT_DIRECTORY)/html
export DOXYGEN_STRIP_PATH := $(DOXYGEN_SOURCE_DIRECTORY)
export DOXYGEN_INCLUDE_PATH := $(shell find $(DOXYGEN_STRIP_PATH) -type d -name dictpch -prune -o -type d -name inc)

LOGO_PATH = $(HOME_WAT)/tools/cwb/www/logo/cluster_site_logo_modern.png
ifeq ($(SITE_CLUSTER), CIT)
LOGO_PATH = $(HOME_WAT)/tools/cwb/www/logo/ligo_virgo_logo_modern.png
endif
ifeq ($(SITE_CLUSTER), ATLAS)
LOGO_PATH = $(HOME_WAT)/tools/cwb/www/logo/atlas_logo_modern.png
endif
ifeq ($(SITE_CLUSTER), VIRTUALBOX)
LOGO_PATH = $(HOME_WAT)/tools/cwb/www/logo/cluster_site_logo_modern.png
endif

WDIR = ${PWD}

define MkDir
	+@[ -d $1 ] || mkdir -p $1
endef

all: filter folders images doxygen

filter:
	`root-config --cxx` -o filter filter.cxx -std=c++11

images:
	$(call MkDir,$(DOXYGEN_IMAGE_PATH))
	cp images/* $(DOXYGEN_IMAGE_PATH)
	cp $(HOME_WAT)/html/cwb-lbanner-modern.png $(DOXYGEN_IMAGE_PATH)
	cp $(LOGO_PATH) $(DOXYGEN_IMAGE_PATH)/logo_right_side.png

doxygen: images filter
	$(call MkDir,$(DOXYGEN_OUTPUT_DIRECTORY))
	$(call MkDir,$(DOXYGEN_EXAMPLE_PATH))
	./makehtmlfooter.sh > htmlfooter.html
	doxygen
	./listlibs.sh
	rm -f listofclass.sh
#	cp ../_build/html/tutorials.html $(DOXYGEN_OUTPUT_DIRECTORY)/html/group__Tutorials.html 

clean:
	rm -rf $(DOXYGEN_OUTPUT_DIRECTORY)/html filter htmlfooter.html 

