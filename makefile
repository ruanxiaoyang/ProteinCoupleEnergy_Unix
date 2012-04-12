all:
	g++ -I. BLOSUM.cpp -o BLOSUM -std=c++0x
	g++ -I. GloAlign.cpp -l pthread -o GloAlign -std=c++0x
	g++ -I. RBLAST.cpp -l pthread -o RBLAST -std=c++0x
	g++ -I. MSA.cpp -l pthread -o MSA -std=c++0x
	g++ -I. CoupleEnergy.cpp -l pthread -o CoupleEnergy -std=c++0x
