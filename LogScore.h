
darray<int> getDayhoff250()
{
    darray<int> D250(20,20,0);
//	    A             R              N             D             C              Q             E             G             H              I              L               K               M               F               P              S              T              W               Y               V 
/*A*/  D250(0,0)=2;  D250(0,1)=-2;  D250(0,2)=0;  D250(0,3)=0;  D250(0,4)=-2;  D250(0,5)=0;  D250(0,6)=0;  D250(0,7)=1;  D250(0,8)=-1;  D250(0,9)=-1;  D250(0,10)=-2;  D250(0,11)=-1;  D250(0,12)=-1;  D250(0,13)=-4;  D250(0,14)=1;  D250(0,15)=1;  D250(0,16)=1;  D250(0,17)=-6;  D250(0,18)=-3;  D250(0,19)=0;
/*R*/  D250(1,0)=-2; D250(1,1)=6;   D250(1,2)=0;  D250(1,3)=-1; D250(1,4)=-4;  D250(1,5)=1;  D250(1,6)=-1; D250(1,7)=-3; D250(1,8)=2;   D250(1,9)=-2;  D250(1,10)=-3;  D250(1,11)=3;   D250(1,12)=0;   D250(1,13)=-4;  D250(1,14)=0;  D250(1,15)=0;  D250(1,16)=-1; D250(1,17)=2;   D250(1,18)=-4;  D250(1,19)=-2;
/*N*/  D250(2,0)=0;  D250(2,1)=0;   D250(2,2)=2;  D250(2,3)=2;  D250(2,4)=-4;  D250(2,5)=1;  D250(2,6)=1;  D250(2,7)=0;  D250(2,8)=2;   D250(2,9)=-2;  D250(2,10)=-3;  D250(2,11)=1;   D250(2,12)=-2;  D250(2,13)=-4;  D250(2,14)=-1; D250(2,15)=1;  D250(2,16)=0;  D250(2,17)=-4;  D250(2,18)=-2;  D250(2,19)=-2;
/*D*/  D250(3,0)=0;  D250(3,1)=-1;  D250(3,2)=2;  D250(3,3)=4;  D250(3,4)=-5;  D250(3,5)=2;  D250(3,6)=3;  D250(3,7)=1;  D250(3,8)=1;   D250(3,9)=-2;  D250(3,10)=-4;  D250(3,11)=0;   D250(3,12)=-3;  D250(3,13)=-6;  D250(3,14)=-1; D250(3,15)=0;  D250(3,16)=0;  D250(3,17)=-7;  D250(3,18)=-4;  D250(3,19)=-2;
/*C*/  D250(4,0)=-2; D250(4,1)=-4;  D250(4,2)=-4; D250(4,3)=-5; D250(4,4)=12;  D250(4,5)=-5; D250(4,6)=-5; D250(4,7)=-3; D250(4,8)=-3;  D250(4,9)=-2;  D250(4,10)=-6;  D250(4,11)=-5;  D250(4,12)=-5;  D250(4,13)=-4;  D250(4,14)=-3; D250(4,15)=0;  D250(4,16)=-2; D250(4,17)=-8;  D250(4,18)=0;   D250(4,19)=-2;
/*Q*/  D250(5,0)=0;  D250(5,1)=1;   D250(5,2)=1;  D250(5,3)=2;  D250(5,4)=-5;  D250(5,5)=4;  D250(5,6)=2;  D250(5,7)=-1; D250(5,8)=3;   D250(5,9)=-2;  D250(5,10)=-2;  D250(5,11)=1;   D250(5,12)=-1;  D250(5,13)=-5;  D250(5,14)=0;  D250(5,15)=-1; D250(5,16)=-1; D250(5,17)=-5;  D250(5,18)=-4;  D250(5,19)=-2;
/*E*/  D250(6,0)=0;  D250(6,1)=-1;  D250(6,2)=1;  D250(6,3)=3;  D250(6,4)=-5;  D250(6,5)=2;  D250(6,6)=4;  D250(6,7)=0;  D250(6,8)=1;   D250(6,9)=-2;  D250(6,10)=-3;  D250(6,11)=0;   D250(6,12)=-2;  D250(6,13)=-5;  D250(6,14)=-1; D250(6,15)=0;  D250(6,16)=0;  D250(6,17)=-7;  D250(6,18)=-4;  D250(6,19)=-2;
/*G*/  D250(7,0)=1;  D250(7,1)=-3;  D250(7,2)=0;  D250(7,3)=1;  D250(7,4)=-3;  D250(7,5)=-1; D250(7,6)=0;  D250(7,7)=5;  D250(7,8)=-2;  D250(7,9)=-3;  D250(7,10)=-4;  D250(7,11)=-2;  D250(7,12)=-3;  D250(7,13)=-5;  D250(7,14)=-1; D250(7,15)=1;  D250(7,16)=0;  D250(7,17)=-7;  D250(7,18)=-5;  D250(7,19)=-1;
/*H*/  D250(8,0)=-1; D250(8,1)=2;   D250(8,2)=2;  D250(8,3)=1;  D250(8,4)=-3;  D250(8,5)=3;  D250(8,6)=1;  D250(8,7)=-2; D250(8,8)=6;   D250(8,9)=-2;  D250(8,10)=-2;  D250(8,11)=0;   D250(8,12)=-2;  D250(8,13)=-2;  D250(8,14)=0;  D250(8,15)=-1; D250(8,16)=-1; D250(8,17)=-3;  D250(8,18)=0;   D250(8,19)=-2;
/*I*/  D250(9,0)=-1; D250(9,1)=-2;  D250(9,2)=-2; D250(9,3)=-25; D250(9,4)=-2;  D250(9,5)=-2; D250(9,6)=-2; D250(9,7)=-3; D250(9,8)=-2;  D250(9,9)=5;   D250(9,10)=2;   D250(9,11)=-2;  D250(9,12)=2;   D250(9,13)=1;   D250(9,14)=-2; D250(9,15)=-1; D250(9,16)=0;  D250(9,17)=-5;  D250(9,18)=-1;  D250(9,19)=4;
/*L*/  D250(10,0)=-2; D250(10,1)=-3;  D250(10,2)=-3; D250(10,3)=-4; D250(10,4)=-6;  D250(10,5)=-2; D250(10,6)=-3; D250(10,7)=-4; D250(10,8)=-2;  D250(10,9)=2;   D250(10,10)=6;   D250(10,11)=-3;  D250(10,12)=4;   D250(10,13)=2;   D250(10,14)=-3; D250(10,15)=-3; D250(10,16)=-2; D250(10,17)=-2;  D250(10,18)=-1;  D250(10,19)=2;
/*K*/  D250(11,0)=-1; D250(11,1)=3;   D250(11,2)=1;  D250(11,3)=0;  D250(11,4)=-5;  D250(11,5)=1;  D250(11,6)=0;  D250(11,7)=-2; D250(11,8)=0;   D250(11,9)=-2;  D250(11,10)=-3;  D250(11,11)=5;   D250(11,12)=0;   D250(11,13)=-5;  D250(11,14)=-1; D250(11,15)=0;  D250(11,16)=0;  D250(11,17)=-3;  D250(11,18)=-4;  D250(11,19)=-2;
/*M*/  D250(12,0)=-1; D250(12,1)=0;   D250(12,2)=-2; D250(12,3)=-3; D250(12,4)=-5;  D250(12,5)=-1; D250(12,6)=-2; D250(12,7)=-3; D250(12,8)=-2;  D250(12,9)=2;   D250(12,10)=4;   D250(12,11)=0;   D250(12,12)=6;   D250(12,13)=0;   D250(12,14)=-2; D250(12,15)=-2; D250(12,16)=-1; D250(12,17)=-4;  D250(12,18)=-2;  D250(12,19)=2;
/*F*/  D250(13,0)=-4; D250(13,1)=-4;  D250(13,2)=-4; D250(13,3)=-6; D250(13,4)=-4;  D250(13,5)=-5; D250(13,6)=-5; D250(13,7)=-5; D250(13,8)=-2;  D250(13,9)=1;   D250(13,10)=2;   D250(13,11)=-5;  D250(13,12)=0;   D250(13,13)=9;   D250(13,14)=-5; D250(13,15)=-3; D250(13,16)=-3; D250(13,17)=0;   D250(13,18)=7;   D250(13,19)=-1;
/*P*/  D250(14,0)=1;  D250(14,1)=0;   D250(14,2)=-1; D250(14,3)=-1; D250(14,4)=-3;  D250(14,5)=0;  D250(14,6)=-1; D250(14,7)=-1; D250(14,8)=0;   D250(14,9)=-2;  D250(14,10)=-3;  D250(14,11)=-1;  D250(14,12)=-2;  D250(14,13)=-5;  D250(14,14)=6;  D250(14,15)=1;  D250(14,16)=0;  D250(14,17)=-6;  D250(14,18)=-5;  D250(14,19)=-1;
/*S*/  D250(15,0)=1;  D250(15,1)=0;   D250(15,2)=1;  D250(15,3)=0;  D250(15,4)=0;   D250(15,5)=-1; D250(15,6)=0;  D250(15,7)=1;  D250(15,8)=-1;  D250(15,9)=-1;  D250(15,10)=-3;  D250(15,11)=0;   D250(15,12)=-2;  D250(15,13)=-3;  D250(15,14)=1;  D250(15,15)=2;  D250(15,16)=1;  D250(15,17)=-2;  D250(15,18)=-3;  D250(15,19)=-1;
/*T*/  D250(16,0)=1;  D250(16,1)=-1;  D250(16,2)=0;  D250(16,3)=0;  D250(16,4)=-2;  D250(16,5)=-1; D250(16,6)=0;  D250(16,7)=0;  D250(16,8)=-1;  D250(16,9)=0;   D250(16,10)=-2;  D250(16,11)=0;   D250(16,12)=-1;  D250(16,13)=-3;  D250(16,14)=0;  D250(16,15)=1;  D250(16,16)=3;  D250(16,17)=-5;  D250(16,18)=-3;  D250(16,19)=0;
/*W*/  D250(17,0)=-6; D250(17,1)=2;   D250(17,2)=-4; D250(17,3)=-7; D250(17,4)=-8;  D250(17,5)=-5; D250(17,6)=-7; D250(17,7)=-7; D250(17,8)=-3;  D250(17,9)=-5;  D250(17,10)=-2;  D250(17,11)=-3;  D250(17,12)=-4;  D250(17,13)=0;   D250(17,14)=-6; D250(17,15)=-2; D250(17,16)=-5; D250(17,17)=17;  D250(17,18)=0;   D250(17,19)=-6;
/*Y*/  D250(18,0)=-3; D250(18,1)=-4;  D250(18,2)=-2; D250(18,3)=-4; D250(18,4)=0;   D250(18,5)=-4; D250(18,6)=-4; D250(18,7)=-5; D250(18,8)=0;   D250(18,9)=-1;  D250(18,10)=-1;  D250(18,11)=-4;  D250(18,12)=-2;  D250(18,13)=7;   D250(18,14)=-5; D250(18,15)=-3; D250(18,16)=-3; D250(18,17)=0;   D250(18,18)=10;  D250(18,19)=-2;
/*V*/  D250(19,0)=0;  D250(19,1)=-2;  D250(19,2)=-2; D250(19,3)=-2; D250(19,4)=-2;  D250(19,5)=-2; D250(19,6)=-2; D250(19,7)=-1; D250(19,8)=-2;  D250(19,9)=4;   D250(19,10)=2;   D250(19,11)=-2;  D250(19,12)=2;   D250(19,13)=-1;  D250(19,14)=-1; D250(19,15)=-1; D250(19,16)=0;  D250(19,17)=-6;  D250(19,18)=-2;  D250(19,19)=4;
   return D250;
}

darray<double> GonnetScorematx()
{
      darray<double> G(20,20,0.0);
//     A                R                N                D                C                Q                E                G                H                I                L                 K                 M                 F                 P                 S                 T                 W                 Y                 V 
/*A*/G(0,0)=2.4;    G(0,1)=-0.6;   G(0,2)=-0.3;   G(0,3)=-0.3;   G(0,4)=0.5;    G(0,5)=-0.2;   G(0,6)=0.0;    G(0,7)=0.5;    G(0,8)=-0.9;   G(0,9)=-0.8;   G(0,10)=-1.2;   G(0,11)=-0.4;   G(0,12)=-0.7;   G(0,13)=-2.3;   G(0,14)=0.3;    G(0,15)=1.1;    G(0,16)=0.6;    G(0,17)=-3.6;   G(0,18)=-2.2;   G(0,19)=0.1;
/*R*/G(1,0)=-0.6;   G(1,1)=4.7;    G(1,2)=0.3;    G(1,3)=-0.3;   G(1,4)=-2.2;   G(1,5)=1.5;    G(1,6)=0.4;    G(1,7)=-1.0;   G(1,8)=0.6;    G(1,9)=-2.4;   G(1,10)=-2.2;   G(1,11)=2.7;    G(1,12)=-1.7;   G(1,13)=-3.2;   G(1,14)=-0.9;   G(1,15)=-0.2;   G(1,16)=-0.2;   G(1,17)=-1.6;   G(1,18)=-1.8;   G(1,19)=-2.0;
/*N*/G(2,0)=-0.3;   G(2,1)=0.3;    G(2,2)=3.8;    G(2,3)=2.2;    G(2,4)=-1.8;   G(2,5)=0.7;    G(2,6)=0.9;    G(2,7)=0.4;    G(2,8)=1.2;    G(2,9)=-2.8;   G(2,10)=-3.0;   G(2,11)=0.8;    G(2,12)=-2.2;   G(2,13)=-3.1;   G(2,14)=-0.9;   G(2,15)=0.9;    G(2,16)=0.5;    G(2,17)=-3.6;   G(2,18)=-1.4;   G(2,19)=-2.2;
/*D*/G(3,0)=-0.3;   G(3,1)=-0.3;   G(3,2)=2.2;    G(3,3)=4.7;    G(3,4)=-3.2;   G(3,5)=0.9;    G(3,6)=2.7;    G(3,7)=0.1;    G(3,8)=0.4;    G(3,9)=-3.8;   G(3,10)=-4.0;   G(3,11)=0.5;    G(3,12)=-3.0;   G(3,13)=-4.5;   G(3,14)=-0.7;   G(3,15)=0.5;    G(3,16)=0.0;    G(3,17)=-5.2;   G(3,18)=-2.8;   G(3,19)=-2.9;
/*C*/G(4,0)=0.5;    G(4,1)=-2.2;   G(4,2)=-1.8;   G(4,3)=-3.2;   G(4,4)=11.5;   G(4,5)=-2.4;   G(4,6)=-3.0;   G(4,7)=-2.0;   G(4,8)=-1.3;   G(4,9)=-1.1;   G(4,10)=-1.5;   G(4,11)=-2.8;   G(4,12)=-0.9;   G(4,13)=-0.8;   G(4,14)=-3.1;   G(4,15)=0.1;    G(4,16)=-0.5;   G(4,17)=-1.0;   G(4,18)=-0.5;   G(4,19)=0.0;
/*Q*/G(5,0)=-0.2;   G(5,1)=1.5;    G(5,2)=0.7;    G(5,3)=0.9;    G(5,4)=-2.4;   G(5,5)=2.7;    G(5,6)=1.7;    G(5,7)=-1.0;   G(5,8)=1.2;    G(5,9)=-1.9;   G(5,10)=-1.6;   G(5,11)=1.5;    G(5,12)=-1.0;   G(5,13)=-2.6;   G(5,14)=-0.2;   G(5,15)=0.2;    G(5,16)=0.0;    G(5,17)=-2.7;   G(5,18)=-1.7;   G(5,19)=-1.5;
/*E*/G(6,0)=0.0;    G(6,1)=0.4;    G(6,2)=0.9;    G(6,3)=2.7;    G(6,4)=-3.0;   G(6,5)=1.7;    G(6,6)=3.6;    G(6,7)=-0.8;   G(6,8)=0.4;    G(6,9)=-2.7;   G(6,10)=-2.8;   G(6,11)=1.2;    G(6,12)=-2.0;   G(6,13)=-3.9;   G(6,14)=-0.5;   G(6,15)=0.2;    G(6,16)=-0.1;   G(6,17)=-4.3;   G(6,18)=-2.7;   G(6,19)=-1.9;
/*G*/G(7,0)=0.5;    G(7,1)=-1.0;   G(7,2)=0.4;    G(7,3)=0.1;    G(7,4)=-2.0;   G(7,5)=-1.0;   G(7,6)=-0.8;   G(7,7)=6.6;    G(7,8)=-1.4;   G(7,9)=-4.5;   G(7,10)=-4.4;   G(7,11)=-1.1;   G(7,12)=-3.5;   G(7,13)=-5.2;   G(7,14)=-1.6;   G(7,15)=0.4;    G(7,16)=-1.1;   G(7,17)=-4.0;   G(7,18)=-4.0;   G(7,19)=-3.3;
/*H*/G(8,0)=-0.9;   G(8,1)=0.6;    G(8,2)=1.2;    G(8,3)=0.4;    G(8,4)=-1.3;   G(8,5)=1.2;    G(8,6)=0.4;    G(8,7)=-1.4;   G(8,8)=6.0;    G(8,9)=-2.2;   G(8,10)=-1.9;   G(8,11)=0.6;    G(8,12)=-1.3;   G(8,13)=-0.1;   G(8,14)=-1.1;   G(8,15)=-0.2;   G(8,16)=-0.3;   G(8,17)=-0.8;   G(8,18)=2.2;    G(8,19)=-2.0;
/*I*/G(9,0)=-0.8;   G(9,1)=-2.4;   G(9,2)=-2.8;   G(9,3)=-3.8;   G(9,4)=-1.1;   G(9,5)=-1.9;   G(9,6)=-2.7;   G(9,7)=-4.5;   G(9,8)=-2.2;   G(9,9)=4.0;    G(9,10)=2.8;    G(9,11)=-2.1;   G(9,12)=2.5;    G(9,13)=1.0;    G(9,14)=-2.6;   G(9,15)=-1.8;   G(9,16)=-0.6;   G(9,17)=-1.8;   G(9,18)=-0.7;   G(9,19)=3.1;
/*L*/G(10,0)=-1.2;  G(10,1)=-2.2;  G(10,2)=-3.0;  G(10,3)=-4.0;  G(10,4)=-1.5;  G(10,5)=-1.6;  G(10,6)=-2.8;  G(10,7)=-4.4;  G(10,8)=-1.9;  G(10,9)=2.8;   G(10,10)=4.0;   G(10,11)=-2.1;  G(10,12)=2.8;   G(10,13)=2.0;   G(10,14)=-2.3;  G(10,15)=-2.1;  G(10,16)=-1.3;  G(10,17)=-0.7;  G(10,18)=0.0;   G(10,19)=1.8;
/*K*/G(11,0)=-0.4;  G(11,1)=2.7;   G(11,2)=0.8;   G(11,3)=0.5;   G(11,4)=-2.8;  G(11,5)=1.5;   G(11,6)=1.2;   G(11,7)=-1.1;  G(11,8)=0.6;   G(11,9)=-2.1;  G(11,10)=-2.1;  G(11,11)=3.2;   G(11,12)=-1.4;  G(11,13)=-3.3;  G(11,14)=-0.6;  G(11,15)=0.1;   G(11,16)=0.1;   G(11,17)=-3.5;  G(11,18)=-2.1;  G(11,19)=-1.7;
/*M*/G(12,0)=-0.7;  G(12,1)=-1.7;  G(12,2)=-2.2;  G(12,3)=-3.0;  G(12,4)=-0.9;  G(12,5)=-1.0;  G(12,6)=-2.0;  G(12,7)=-3.5;  G(12,8)=-1.3;  G(12,9)=2.5;   G(12,10)=2.8;   G(12,11)=-1.4;  G(12,12)=4.3;   G(12,13)=1.6;   G(12,14)=-2.4;  G(12,15)=-1.4;  G(12,16)=-0.6;  G(12,17)=-1.0;  G(12,18)=-0.2;  G(12,19)=1.6;
/*F*/G(13,0)=-2.3;  G(13,1)=-3.2;  G(13,2)=-3.1;  G(13,3)=-4.5;  G(13,4)=-0.8;  G(13,5)=-2.6;  G(13,6)=-3.9;  G(13,7)=-5.2;  G(13,8)=-0.1;  G(13,9)=1.0;   G(13,10)=2.0;   G(13,11)=-3.3;  G(13,12)=1.6;   G(13,13)=7.0;   G(13,14)=-3.8;  G(13,15)=-2.8;  G(13,16)=-2.2;  G(13,17)=3.6;   G(13,18)=5.1;   G(13,19)=0.1;
/*P*/G(14,0)=0.3;   G(14,1)=-0.9;  G(14,2)=-0.9;  G(14,3)=-0.7;  G(14,4)=-3.1;  G(14,5)=-0.2;  G(14,6)=-0.5;  G(14,7)=-1.6;  G(14,8)=-1.1;  G(14,9)=-2.6;  G(14,10)=-2.3;  G(14,11)=-0.6;  G(14,12)=-2.4;  G(14,13)=-3.8;  G(14,14)=7.6;   G(14,15)=0.4;   G(14,16)=0.1;   G(14,17)=-5.0;  G(14,18)=-3.1;  G(14,19)=-1.8;
/*S*/G(15,0)=1.1;   G(15,1)=-0.2;  G(15,2)=0.9;   G(15,3)=0.5;   G(15,4)=0.1;   G(15,5)=0.2;   G(15,6)=0.2;   G(15,7)=0.4;   G(15,8)=-0.2;  G(15,9)=-1.8;  G(15,10)=-2.1;  G(15,11)=0.1;   G(15,12)=-1.4;  G(15,13)=-2.8;  G(15,14)=0.4;   G(15,15)=2.2;   G(15,16)=1.5;   G(15,17)=-3.3;  G(15,18)=-1.9;  G(15,19)=-1.0;
/*T*/G(16,0)=0.6;   G(16,1)=-0.2;  G(16,2)=0.5;   G(16,3)=0.0;   G(16,4)=-0.5;  G(16,5)=0.0;   G(16,6)=-0.1;  G(16,7)=-1.1;  G(16,8)=-0.3;  G(16,9)=-0.6;  G(16,10)=-1.3;  G(16,11)=0.1;   G(16,12)=-0.6;  G(16,13)=-2.2;  G(16,14)=0.1;   G(16,15)=1.5;   G(16,16)=2.5;   G(16,17)=-1.8;  G(16,18)=-1.9;  G(16,19)=0.0;
/*W*/G(17,0)=-3.6;  G(17,1)=-1.6;  G(17,2)=-3.6;  G(17,3)=-5.2;  G(17,4)=-1.0;  G(17,5)=-2.7;  G(17,6)=-4.3;  G(17,7)=-4.0;  G(17,8)=-0.8;  G(17,9)=-1.8;  G(17,10)=-0.7;  G(17,11)=-3.5;  G(17,12)=-1.0;  G(17,13)=3.6;   G(17,14)=-5.0;  G(17,15)=-3.3;  G(17,16)=-1.8;  G(17,17)=14.2;  G(17,18)=4.1;   G(17,19)=-2.6;
/*Y*/G(18,0)=-2.2;  G(18,1)=-1.8;  G(18,2)=-1.4;  G(18,3)=-2.8;  G(18,4)=-0.5;  G(18,5)=-1.7;  G(18,6)=-2.7;  G(18,7)=-4.0;  G(18,8)=2.2;   G(18,9)=-0.7;  G(18,10)=0.0;   G(18,11)=-2.1;  G(18,12)=-0.2;  G(18,13)=5.1;   G(18,14)=-3.1;  G(18,15)=-1.9;  G(18,16)=-1.9;  G(18,17)=4.1;   G(18,18)=7.8;   G(18,19)=-1.1;
/*V*/G(19,0)=0.1;   G(19,1)=-2.0;  G(19,2)=-2.2;  G(19,3)=-2.9;  G(19,4)=0.0;   G(19,5)=-1.5;  G(19,6)=-1.9;  G(19,7)=-3.3;  G(19,8)=-2.0;  G(19,9)=3.1;   G(19,10)=1.8;   G(19,11)=-1.7;  G(19,12)=1.6;   G(19,13)=0.1;   G(19,14)=-1.8;  G(19,15)=-1.0;  G(19,16)=0.0;   G(19,17)=-2.6;  G(19,18)=-1.1;  G(19,19)=3.4;
    return G;
}


int aaid(char & _aa)
{
	if((int)_aa<=70 || ((int)_aa>=97 && (int)_aa<=102) )
	{
		if(_aa=='A' || _aa=='a')
			return 0;
		if(_aa=='C' || _aa=='c')
			return 4;
		if(_aa=='D' || _aa=='d')
			return 3;
		if(_aa=='E' || _aa=='e')
			return 6;
		if(_aa=='F' || _aa=='f')
			return 13;
		if(_aa=='-')
			return -4;
		if(_aa=='B' || _aa=='b')
			return 20;
	}
	else if((int)_aa<=76 || ((int)_aa>=103 && (int)_aa<=108))
	{
		if(_aa=='G' || _aa=='g')
			return 7;
		if(_aa=='H' || _aa=='h')
			return 8;
		if(_aa=='I' || _aa=='i')
			return 9;
		if(_aa=='K' || _aa=='k')
			return 11;
		if(_aa=='L' || _aa=='l')
			return 10;
	}
	else if((int)_aa<=82 || ((int)_aa>=109 && (int)_aa<=114))
	{
		if(_aa=='M' || _aa=='m')
			return 12;
		if(_aa=='N' || _aa=='n')
			return 2;
		if(_aa=='P' || _aa=='p')
			return 14;
		if(_aa=='Q' || _aa=='q')
			return 5;
		if(_aa=='R' || _aa=='r')
			return 1;
	}
	else
	{
		if(_aa=='S' || _aa=='s')
			return 15;
		if(_aa=='T' || _aa=='t')
			return 16;
		if(_aa=='V' || _aa=='v')
			return 19;
		if(_aa=='W' || _aa=='w')
			return 17;
		if(_aa=='Y' || _aa=='y')
			return 18;
		if(_aa=='Z' || _aa=='z')
			return 21;
		if(_aa=='X' || _aa=='x')
			return 22;
	}
	return rand()%20;
}
char idaa(const int & id)
{
	if(id<=4)
	{
		if(id==0) return 'A';
		if(id==1) return 'R';
		if(id==2) return 'N';
		if(id==3) return 'D';
		if(id==4) return 'C';
		if(id==-4) return '-';
	}
	else if(id<=9)
	{
		if(id==5) return 'Q';
		if(id==6) return 'E';
		if(id==7) return 'G';
		if(id==8) return 'H';
		if(id==9) return 'I';
	}
	else if(id<=14)
	{
		if(id==10) return 'L';
		if(id==11) return 'K';
		if(id==12) return 'M';
		if(id==13) return 'F';
		if(id==14) return 'P';
	}
	else if(id<=19)
	{
		if(id==15) return 'S';
		if(id==16) return 'T';
		if(id==17) return 'W';
		if(id==18) return 'Y';
		if(id==19) return 'V';
	}
	else
	{
		if(id==20) return 'B';
		if(id==21) return 'Z';
		if(id==22) return 'X';
	}
	return 'X';
}
template <class type>
type SCORE(darray<type> _scorematx,char & _a,char & _b)
{
	int i=aaid(_a);
	int j=aaid(_b);
	if(i>=0 && j>=0)
		return _scorematx(i,j);
	else if(i==-3 && j>=0)
        return _scorematx(j,j);
	else if(i>=0 && j==-3)
		return _scorematx(i,i);
	else if((i==-1 && j==-1) || (i==-1 && j==3) || (i==3 && j==-1) || (i==-3 && j==-1) || (i==-1 && j==-3) )
		return _scorematx(3,3);
	else if((i==-1 && j==2) || (i==2 && j==-1))
		return _scorematx(2,2);
	else if((i==-2 && j==-2) || (i==-2 && j==6) || (i==6 && j==-2) || (i==-3 && j==-2) || (i==-2 && j==-3) )
		return _scorematx(6,6);
	else if((i==-2 && j==5) || (i==5 && j==-2))
		return _scorematx(5,5);
	else
		return -10;
}

sarray<char> aa20()
{
	sarray<char> aas(20);
	aas[0]='A';aas[1]='R';aas[2]='N';aas[3]='D';aas[4]='C';
	aas[5]='Q';aas[6]='E';aas[7]='G';aas[8]='H';aas[9]='I';
	aas[10]='L';aas[11]='K';aas[12]='M';aas[13]='F';aas[14]='P';
	aas[15]='S';aas[16]='T';aas[17]='W';aas[18]='Y';aas[19]='V';
	return aas;
}
