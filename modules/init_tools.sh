git clone https://github.com/suhrig/arriba
cd ./arriba && make
cd ..
git clone https://github.com/alexdobin/STAR
cd ./STAR/source && make
cd ../../
git clone https://github.com/DavidsonGroup/flexiplex
cd ./flexiplex && make
cd ../barcodes/ && gunzip 3M-february-2018.txt.gz
