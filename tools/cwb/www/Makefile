all : 	modern 

classic : clean xclassic index
modern :  clean xmodern index

xclassic :  
	ln -sf cedindex_header_classic.html cedindex_header.html 

xmodern : 
	ln -sf cedindex_header_modern.html cedindex_header.html 

index : 
	if [ ! -d index ]; then mkdir index ; fi;
	cp ced.css index/. 
	cp ced.js index/. 
	cp logo/* index/. 
	root -l -b LoopCreateIndexCED.C 
	rm -f cedindex_header.html

clean : 
	rm -rf index
	rm -f cedindex_header.html

