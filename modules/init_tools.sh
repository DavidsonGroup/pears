git clone https://github.com/suhrig/arriba
cd ./arriba && make
cd ..
git clone https://github.com/alexdobin/STAR
cd ./STAR/source && make
cd ../../
git clone https://github.com/DavidsonGroup/flexiplex
cd ./flexiplex && make
cd ../barcodes/ && gunzip 3M-february-2018.txt.gz
cd ..
wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2
tar -xf samtools-1.18.tar.bz2
cd samtools-1.18.tar.bz2
make
