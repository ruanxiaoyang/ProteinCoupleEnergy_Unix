all:
	g++ -I. BLOSUM.cpp -o BLOSUM -std=c++0x
	g++ -I. GloAlign.cpp -l pthread -o GloAlign -std=c++0x
	g++ -I. RBLAST.cpp -l pthread -o RBLAST -std=c++0x
	g++ -I. MSA.cpp -l pthread -o MSA -std=c++0x
	g++ -I. CoupleEnerge.cpp -l pthread -o CoupleEnerge -std=c++0x
