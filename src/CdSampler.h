/*
 * CdSampler.h
 *
 *  Created on: Jan 2, 2014
 *      Author: ianfellows
 */

#ifndef CDSAMPLER_H_
#define CDSAMPLER_H_

#include "Rcpp.h"
#include <cmath>
#include <vector>
#include <assert.h>
#include <memory>
#include <boost/shared_ptr.hpp>
#include <boost/unordered_set.hpp>
#include "MetropolisHastings.h"
#include "DyadToggle.h"
#include "VertexToggle.h"
#include "VertexToggles.h"
#include "ToggleController.h"
#include "ShallowCopyable.h"

namespace ernm {


template<class Engine>
class DefaultCd{

protected:
	typedef typename BinaryNet<Engine>::NeighborIterator NeighborIterator;
	typedef typename boost::unordered_set< std::pair<int,int>, PairHash >  DyadSet;
	//typedef typename std::tr1::unordered_set< int >  NbrSet;
	typedef std::vector< int >  NbrSet;
	typedef boost::shared_ptr< AbstractDyadToggle<Engine> > DyadTogglePtr;
	typedef boost::shared_ptr< AbstractVertexToggle<Engine> > VertexTogglePtr;
	boost::shared_ptr< BinaryNet<Engine> > net;

	//Neighborhood<Engine> nt;
	NodeTieDyad<Engine> ntd;
	//RandomDyad<Engine> ntd;

	std::vector< std::pair<int,int> > toggle;

	double p;
	double alpha;
	double n;
	int activeNode;
	int lastMethod;

	std::vector<int> order;
	std::set<int> from;
	std::set<int> to;
	Set nbr;
	int n1;
	int n2;

	int ss;
	DyadSet dyads;
	NbrSet nbrs;
	//std::vector<int> nbrs;
	double lr;

	void addToggle(std::pair<int,int> tog){
		if(tog.first < tog.second)
			dyads.insert(tog);
		else
			dyads.insert(std::pair<int,int>(tog.second,tog.first));
	}

public:


	DefaultCd(){
		n=0;
		alpha = .5;
		p = 1.0 / 3.0;
		activeNode = -1;
		lr = 0.0;
		lastMethod = 0;
		n1 = n2 = -1;
		ss=20;
	}

	DefaultCd(Rcpp::List l){
		n=0.0;
		alpha = .5;
		p = 1.0 / 3.0;
		activeNode = -1;
		lr = 0.0;
		lastMethod = 0;
		n1 = n2 = -1;
		try{
			ss = as< int >(l(0));
		}catch(...){
			ss = 20;
		}
		std::cout << ss;
	}

	DefaultCd( BinaryNet<Engine> & network){
		boost::shared_ptr< BinaryNet<Engine> > n(new  BinaryNet<Engine> (network));
		setNetwork(n);
		this->n = 0.0;
		alpha = .5;
		p = .8;
		activeNode = -1;
		lr = 0.0;
		lastMethod = 0;
	}

	DefaultCd( BinaryNet<Engine> & network, int sampSize){
		boost::shared_ptr< BinaryNet<Engine> > n(new  BinaryNet<Engine> (network));
		setNetwork(n);
		this->n = 0.0;
		alpha = .5;
		p = .8;
		activeNode = -1;
		lr = 0.0;
		lastMethod = 0;
		ss = sampSize;
	}

	virtual ~DefaultCd(){}

	void setOrd(){
		int ord[1270] = //{16, 14, 3, 1, 8, 33, 21, 25, 15, 18, 11, 28, 27, 31, 30, 23,
						//9, 19, 26, 35, 13, 24, 12, 29, 2, 6, 20, 7, 22, 0, 10, 17, 34,
						//32, 5, 4};
		//{13, 24, 2, 6, 20, 0, 10, 7, 22, 12, 29, 17, 34, 32, 4, 5, 23,
		//9, 19, 26, 35, 27, 30, 31, 16, 14, 3, 1, 8, 33, 21, 25, 15, 18,
		//11, 28};
		//{14, 26, 28, 33, 35, 25, 23, 22, 21, 20, 19, 12, 10, 0, 7, 18,
		//		16, 15, 11, 9, 8, 1, 3, 32, 31, 30, 29, 17, 5, 2, 4, 34, 27,
		//		24, 6, 13};
		{2, 9, 14, 16, 17, 19, 26, 31, 38, 40, 43, 47, 48, 49, 61, 64,
				68, 69, 75, 77, 85, 86, 87, 88, 94, 96, 98, 105, 106, 109, 114,
				118, 126, 131, 133, 141, 143, 145, 147, 159, 160, 161, 162, 165,
				173, 175, 176, 178, 182, 185, 186, 190, 191, 192, 202, 204, 209,
				214, 215, 224, 240, 241, 243, 247, 249, 251, 252, 255, 260, 263,
				265, 268, 271, 273, 274, 278, 283, 291, 292, 297, 299, 303, 304,
				306, 308, 310, 311, 316, 320, 321, 323, 325, 329, 331, 335, 339,
				341, 342, 343, 345, 347, 348, 355, 357, 359, 361, 363, 366, 370,
				373, 381, 384, 387, 390, 396, 400, 403, 407, 418, 427, 428, 430,
				432, 434, 437, 440, 444, 446, 448, 450, 451, 452, 455, 456, 458,
				462, 466, 469, 470, 472, 473, 478, 479, 482, 491, 494, 495, 500,
				503, 506, 514, 516, 518, 520, 536, 539, 542, 543, 547, 550, 553,
				554, 558, 568, 585, 586, 587, 589, 591, 602, 605, 611, 616, 621,
				627, 628, 634, 645, 650, 660, 661, 664, 666, 670, 676, 678, 683,
				686, 691, 693, 695, 706, 707, 719, 722, 725, 727, 741, 743, 745,
				747, 748, 750, 751, 755, 756, 757, 763, 765, 770, 773, 774, 775,
				785, 786, 787, 788, 795, 797, 800, 802, 806, 808, 811, 812, 817,
				818, 821, 830, 835, 836, 840, 842, 849, 851, 852, 854, 859, 860,
				861, 866, 871, 872, 878, 881, 891, 892, 894, 897, 898, 901, 905,
				911, 913, 918, 919, 922, 923, 933, 937, 944, 957, 962, 970, 972,
				974, 976, 978, 981, 982, 983, 986, 987, 990, 991, 997, 999, 1005,
				1008, 1011, 1012, 1013, 1021, 1030, 1034, 1035, 1040, 1044, 1048,
				1055, 1056, 1057, 1061, 1063, 1064, 1080, 1082, 1085, 1088, 1089,
				1091, 1093, 1098, 1100, 1102, 1103, 1104, 1105, 1106, 1108, 1109,
				1122, 1125, 1127, 1129, 1131, 1134, 1136, 1138, 1141, 1143, 1145,
				1146, 1151, 1153, 1157, 1161, 1163, 1164, 1167, 1169, 1170, 1172,
				1175, 1176, 1178, 1179, 1184, 1186, 1188, 1192, 1194, 1199, 1200,
				1203, 1213, 1217, 1219, 1222, 1229, 1232, 1235, 1236, 1246, 1248,
				1253, 1255, 1258, 1260, 1266, 1268, 3, 12, 25, 28, 33, 35, 37,
				46, 50, 52, 55, 60, 73, 78, 81, 82, 99, 103, 110, 116, 117, 119,
				128, 139, 144, 146, 148, 149, 150, 151, 155, 157, 158, 163, 166,
				169, 172, 180, 181, 187, 188, 203, 206, 207, 210, 211, 217, 220,
				221, 227, 231, 233, 237, 238, 248, 261, 264, 266, 270, 277, 279,
				280, 282, 284, 288, 290, 295, 298, 300, 312, 319, 324, 326, 333,
				336, 340, 346, 351, 356, 360, 367, 378, 391, 393, 397, 401, 404,
				411, 414, 417, 422, 423, 426, 429, 443, 445, 453, 459, 460, 464,
				467, 471, 480, 483, 485, 486, 488, 492, 499, 508, 510, 511, 513,
				515, 517, 519, 525, 526, 535, 537, 541, 544, 549, 557, 560, 575,
				577, 580, 582, 583, 588, 593, 597, 599, 608, 609, 610, 614, 615,
				625, 630, 632, 633, 640, 644, 648, 649, 651, 671, 674, 688, 689,
				696, 697, 699, 700, 704, 705, 708, 710, 721, 724, 736, 752, 753,
				758, 760, 768, 769, 772, 778, 781, 782, 783, 790, 793, 796, 799,
				803, 807, 809, 813, 815, 822, 824, 827, 828, 831, 833, 834, 837,
				838, 841, 843, 845, 846, 850, 856, 857, 863, 864, 868, 869, 870,
				877, 879, 887, 888, 889, 890, 902, 907, 909, 912, 914, 915, 917,
				921, 925, 934, 935, 943, 945, 946, 947, 954, 956, 958, 960, 961,
				963, 967, 969, 971, 973, 977, 979, 984, 994, 996, 1002, 1003,
				1007, 1014, 1016, 1019, 1024, 1028, 1041, 1042, 1049, 1053, 1058,
				1066, 1069, 1070, 1071, 1077, 1083, 1084, 1087, 1090, 1092, 1094,
				1095, 1096, 1111, 1112, 1114, 1118, 1119, 1123, 1124, 1128, 1130,
				1133, 1135, 1142, 1144, 1147, 1149, 1150, 1154, 1156, 1160, 1177,
				1181, 1182, 1189, 1191, 1202, 1205, 1207, 1212, 1220, 1221, 1223,
				1228, 1230, 1231, 1233, 1234, 1247, 1249, 1252, 1257, 1259, 1262,
				1264, 1, 4, 5, 7, 11, 18, 21, 23, 30, 32, 34, 41, 44, 51, 53,
				54, 56, 58, 59, 66, 67, 72, 76, 79, 80, 83, 84, 91, 101, 104,
				108, 113, 115, 120, 121, 122, 127, 129, 130, 136, 137, 138, 140,
				142, 154, 156, 168, 170, 174, 184, 194, 195, 196, 197, 198, 199,
				201, 208, 216, 219, 222, 223, 226, 232, 234, 235, 239, 242, 244,
				254, 257, 258, 259, 272, 275, 276, 285, 286, 289, 293, 294, 296,
				301, 305, 309, 315, 317, 318, 322, 328, 332, 334, 337, 338, 344,
				352, 358, 362, 365, 368, 369, 375, 379, 380, 382, 385, 392, 394,
				398, 399, 408, 409, 410, 419, 421, 425, 435, 439, 441, 442, 447,
				454, 457, 463, 465, 468, 474, 475, 477, 487, 490, 496, 497, 501,
				502, 505, 512, 521, 527, 528, 529, 532, 533, 534, 545, 546, 548,
				551, 555, 559, 561, 562, 564, 565, 569, 570, 571, 574, 576, 579,
				581, 592, 594, 596, 600, 601, 603, 604, 606, 612, 613, 619, 623,
				624, 626, 629, 635, 637, 638, 639, 641, 643, 647, 653, 659, 662,
				663, 665, 668, 677, 679, 680, 681, 685, 687, 701, 702, 709, 712,
				713, 714, 717, 718, 720, 723, 726, 728, 729, 730, 731, 733, 740,
				749, 762, 764, 767, 780, 789, 791, 792, 798, 804, 805, 810, 814,
				819, 820, 823, 825, 829, 832, 839, 844, 848, 858, 862, 865, 867,
				873, 880, 882, 883, 886, 893, 900, 903, 904, 906, 908, 910, 920,
				924, 927, 928, 929, 931, 936, 938, 939, 940, 941, 948, 949, 950,
				951, 953, 959, 965, 966, 968, 975, 985, 988, 989, 992, 995, 998,
				1000, 1010, 1015, 1018, 1025, 1026, 1031, 1032, 1036, 1038, 1039,
				1045, 1046, 1050, 1051, 1052, 1054, 1062, 1067, 1068, 1072, 1073,
				1074, 1075, 1078, 1079, 1101, 1107, 1113, 1115, 1116, 1120, 1121,
				1126, 1132, 1137, 1148, 1155, 1158, 1162, 1165, 1173, 1174, 1183,
				1187, 1190, 1197, 1198, 1201, 1204, 1210, 1211, 1215, 1224, 1226,
				1227, 1238, 1239, 1240, 1241, 1242, 1244, 1245, 1250, 1256, 1261,
				1263, 1267, 1269, 0, 6, 8, 10, 13, 15, 20, 22, 24, 27, 29, 36,
				39, 42, 45, 57, 62, 63, 65, 70, 71, 74, 89, 90, 92, 93, 95, 97,
				100, 102, 107, 111, 112, 123, 124, 125, 132, 134, 135, 152, 153,
				164, 167, 171, 177, 179, 183, 189, 193, 200, 205, 212, 213, 218,
				225, 228, 229, 230, 236, 245, 246, 250, 253, 256, 262, 267, 269,
				281, 287, 302, 307, 313, 314, 327, 330, 349, 350, 353, 354, 364,
				371, 372, 374, 376, 377, 383, 386, 388, 389, 395, 402, 405, 406,
				412, 413, 415, 416, 420, 424, 431, 433, 436, 438, 449, 461, 476,
				481, 484, 489, 493, 498, 504, 507, 509, 522, 523, 524, 530, 531,
				538, 540, 552, 556, 563, 566, 567, 572, 573, 578, 584, 590, 595,
				598, 607, 617, 618, 620, 622, 631, 636, 642, 646, 652, 654, 655,
				656, 657, 658, 667, 669, 672, 673, 675, 682, 684, 690, 692, 694,
				698, 703, 711, 715, 716, 732, 734, 735, 737, 738, 739, 742, 744,
				746, 754, 759, 761, 766, 771, 776, 777, 779, 784, 794, 801, 816,
				826, 847, 853, 855, 874, 875, 876, 884, 885, 895, 896, 899, 916,
				926, 930, 932, 942, 952, 955, 964, 980, 993, 1001, 1004, 1006,
				1009, 1017, 1020, 1022, 1023, 1027, 1029, 1033, 1037, 1043, 1047,
				1059, 1060, 1065, 1076, 1081, 1086, 1097, 1099, 1110, 1117, 1139,
				1140, 1152, 1159, 1166, 1168, 1171, 1180, 1185, 1193, 1195, 1196,
				1206, 1208, 1209, 1214, 1216, 1218, 1225, 1237, 1243, 1251, 1254,
				1265};
		order = std::vector<int>(ord,ord+net->size());
		//order.clear();
		//for(int i=0;i<net->size();i++)
		//	order.push_back(i);
	}

	void initialize(){
		toggle.resize(1,std::make_pair(-1,-1));
		dyads.clear();
		nbrs.clear();
		from.clear();
		to.clear();
		activeNode = -1;
		lastMethod = 0;
		n = 0.0;
		n1 = n2 = -1;
		order.clear();
		//nt.initialize();
		ntd.initialize();
	}

	void setNetwork(const boost::shared_ptr< BinaryNet<Engine> > n){
		net = n;
		//nt.setNetwork(n);
		ntd.setNetwork(n);
		setOrd();
	}

	void togglesAccepted(bool apply){
		n++;
		if(apply){
			int from = toggle.at(0).first;
			int to = toggle.at(0).second;
			if(net->hasEdge(from,to)){
				nbr.insert(to);
			}else{
				nbr.erase(to);
			}
		}
	}

	inline int randWrapper(const int n) { return floor(unif_rand()*n); }

	inline void generate(){
		if(activeNode<0)
			activeNode = floor(Rf_runif(0.0,net->size()));

		//int ss = floor(max(36.5,1.0 + sqrt(n)));
		//int ss = 2;
		/*if(from.size()==0){
			while(from.size()<ss)
				from.insert(floor(Rf_runif(0.0,net->size())));
			while(to.size()<ss)
				to.insert(floor(Rf_runif(0.0,net->size())));
		}
		int i1 = 0;
		int i2 = 0;
		while(i1 == i2){
			i1 = floor(Rf_runif(0.0,ss));
			i2 = floor(Rf_runif(0.0,ss));
			set<int>::iterator it = from.begin();
			advance(it,i1);
			i1 = *it;
			it = to.begin();
			advance(it,i2);
			i2 = *it;
		}
		toggle.at(0) = std::make_pair(i1,i2);*/

		/*if(n1 < 0){
			n1 = floor(Rf_runif(0.0,net->size() ));
			n2 = floor(Rf_runif(0.0, net->size() - 1));
			if(n2 >= n1)
				n2++;
		}*/
		if(order.size()==0){
			sampleWithoutReplacement(net->size(),ss,order);
			n1 = floor(Rf_runif(0.0,ss));
			nbr.clear();
			for(int i=0;i<order.size();i++){
				if(i==n1)
					continue;
				if(net->hasEdge(order.at(n1),order.at(i))){
					nbr.insert(order.at(i));
				}
			}
		}
		int node = order.at(n1);
		int ndyads = order.size()-1;
		if(ndyads < 0){
			n2 = floor(Rf_runif(0.0,ss-1));
			if(n2>=n1)
				n2++;
			toggle.at(0) = std::make_pair(order.at(n1),order.at(n2));
			lr = 0.0;
			return;
		}
		int nedges = nbr.size();
		bool pickEdge = Rf_runif(0.0,1.0)>0.5;
		pickEdge = false;
		bool hasEdge;
		if(nedges==0)
			pickEdge=false;
		int neighbor;
		if(pickEdge){
			hasEdge=true;
			int nbrIndex = floor(Rf_runif(0.0,(double)nedges));
			Set::iterator it = nbr.begin();
			std::advance(it,nbrIndex);
			neighbor = *it;
			toggle.at(0) = std::make_pair(node,neighbor);
			assert(node != neighbor);

		}else{
			neighbor = floor(Rf_runif(0.0,order.size() - 1.0));
			if(neighbor>=n1)
				neighbor++;
			neighbor = order.at(neighbor);
			toggle.at(0) = std::make_pair(node,neighbor);
			hasEdge = net->hasEdge(node,neighbor);
			assert(node != neighbor);
		}

		double tForward, tReverse;
		if(hasEdge){
			if(nedges<1.5)
				tReverse = 1.0/ndyads;
			else
				tReverse = 0.5/ndyads;
			tForward = 0.5/nedges + 0.5/ndyads;

		}else{
			if(nedges<.5)
				tForward = 1.0/ndyads;
			else
				tForward = 0.5/ndyads;
			tReverse = 0.5/(nedges + 1.0) + 0.5/ndyads;

		}
		lr = log( tReverse/tForward );
		lr = 0.0;
		return;



		if(order.size()==0){
			sampleWithoutReplacement(net->size(),ss,order);
		}
		n1 = floor(Rf_runif(0.0,ss));
		n2 = floor(Rf_runif(0.0,std::max(0.1,ss-1.0)));
		if(n2>=n1)
			n2++;
		n1 = order.at(n1);
		if(ss > 1)
			n2 = order.at(n2);
		if(ss > 1 && Rf_runif(0.0,1.0) < .05){
			toggle.at(0) = std::make_pair(n1,n2);
			lr = 0.0;
		}else{
			if(ss == 1 || Rf_runif(0.0,1.0) < .5){
				ntd.generate(n1);
			}else{
				ntd.generate(n2);
			}

			toggle = ntd.dyadToggles();
			lr = ntd.logRatio();
		}
		return;
		/*int i1 = n1 + floor(Rf_runif(0.0,ss));
		int i2 = n2 + floor(Rf_runif(0.0,ss));
		if(i1>=net->size())
			i1 = i1 - net->size();
		if(i2>=net->size() - 1)
			i2 = i2 - net->size() + 1;
		if(i2 >= i1)
			i2++;
		if(i1==i2)
			Rf_error(std::string("cd toggle").c_str() );
		toggle.at(0) = std::make_pair(order.at(i1),order.at(i2));
		lr = 0.0;*/


/*		NeighborIterator it;
		NeighborIterator itend;
		if(net->isDirected()){
			it = net->outBegin(activeNode);
			itend = net->outEnd(activeNode);
		}else{
			it = net->begin(activeNode);
			itend = net->end(activeNode);
		}
*/
		if(nbrs.size() == 0){
			/*int ss = 4;
			nbrs = std::vector<int>(ss,-1);
			for(int i=0;i<ss;i++){
				nbrs[i] = floor(Rf_runif(0.0,net->size()-1));
				if(nbrs[i]>=activeNode)
					nbrs[i]++;
			}
			std::sort(nbrs.begin(),nbrs.end());
			int i=0;
			while(i<(ss-1)){
				if(nbrs[i] == nbrs[i+1]){
					nbrs[i] = floor(Rf_runif(0.0,net->size()-1));
					if(nbrs[i]>=activeNode)
						nbrs[i]++;
					std::sort(nbrs.begin(),nbrs.end());
					i = 0;
				}else
					i++;
			}*/

			/*copy(net->begin(activeNode),net->end(activeNode),back_inserter(nbrs));
			if(nbrs.size() == 0){
				nbrs.push_back(floor(Rf_runif(0.0,nbrs.size())));
			}
			if(nbrs.size() < 2){
				int ind = floor(Rf_runif(0.0,nbrs.size()-1));
				if(ind >= nbrs.at(0))
					ind++;
				nbrs.push_back(ind);
			}
			*/
		}
		/*int ind = floor(Rf_runif(0.0,nbrs.size()));
		int ind2 = floor(Rf_runif(0.0,nbrs.size()-1));
		if(ind2>=ind)
			ind2++;
		NbrSet::iterator it = nbrs.begin();
		advance(it,ind);
		NbrSet::iterator it2 = nbrs.begin();
		advance(it2,ind2);
		//toggle[0] = std::pair<int,int>(*it,*it2);
		toggle[0] = std::pair<int,int>(*it,activeNode);
		lr = 0.0;
		lastMethod = 1;
		*/
		/*while(dyads.size()<20){
			dyads.insert(net->randomDyad());
		}
		int ind = floor(Rf_runif(0.0,dyads.size()));
		std::tr1::unordered_set< std::pair<int,int>, PairHash>::iterator it = dyads.begin();
		std::advance(it,ind);
		toggle.at(0) = *it;*/
		return;

		if(dyads.size()>1 && n <= alpha * pow(dyads.size(), 2.0)){
			//select a random dyad already sampled.
			int ind = floor(Rf_runif(0.0,dyads.size()));
			boost::unordered_set< std::pair<int,int>, PairHash>::iterator it = dyads.begin();
			std::advance(it,ind);
			toggle[0] = *it;

			//half the time look at the reverse (only affects directed networks)
			if(Rf_runif(0.0,1.0)<.5){
				int tmp = toggle[0].first;
				toggle[0]. first = toggle[0].second;
				toggle[0].second = tmp;
			}
			lr = 0.0;
			lastMethod = 0;
		}else{
			//select a (possibly) new dyad
			if(nbrs.size()>1 && (lastMethod == 0 || lastMethod == 2)){
				//select among nodes from whom a activeNode <-> node link
				//has already been proposed, completing the triangle.
				int ind = floor(Rf_runif(0.0,nbrs.size()));
				int ind2 = floor(Rf_runif(0.0,nbrs.size()-1));
				if(ind2>=ind)
					ind2++;
				NbrSet::iterator it = nbrs.begin();
				advance(it,ind);
				NbrSet::iterator it2 = nbrs.begin();
				advance(it2,ind2);
				toggle[0] = std::pair<int,int>(*it,*it2);
				lr = 0.0;
				lastMethod = 1;

			}else if(lastMethod == 1 || nbrs.size() < 2){
				//select via nodal tie-dyad centered at active node
				//ntd.generate(activeNode);
				//ntd.generate();
				int to = floor(Rf_runif(0.0,net->size() - 1));
				if(to >= activeNode)
					to++;
				toggle.at(0) = std::make_pair(activeNode,to);
				lastMethod = 2;
			}else{
/*				NeighborIterator it;
				int degree;
				if(!net->isDirected()){
					it = net->begin(activeNode);
					degree = net->degree(activeNode);
				}else{
					it = net->outBegin(activeNode);
					degree = net->outdegree(activeNode);
				}
				int ind = floor(Rf_runif(0.0,degree));
				std::advance(it,ind);
				int intNode =*it;

				int to;
				if(net->degree(intNode)<2){
					to = floor(Rf_runif(0.0,net->size()-1));
					if(to>=activeNode)
						to++;
				}else{

				}*/
			}
		}
	}

	inline std::vector< std::pair<int,int> >& dyadToggles(){
		if(lastMethod == 0){
			return toggle;
		}else if(lastMethod == 1){
			//return nt.dyadToggles();
			return toggle;
		}else{
			//return ntd.dyadToggles();
			return toggle;
		}
	}


	inline double logRatio(){
		if(lastMethod == 0){
			return lr;
		}else if(lastMethod == 1){
			//return nt.logRatio();
			return lr;
		}else{
			//return ntd.logRatio();
			return lr;
		}
	}

	inline std::string name(){
		return "DefaultCd";
	}
};

typedef DyadToggle<Directed, DefaultCd<Directed> > DirectedDefaultCdToggle;
typedef DyadToggle<Undirected, DefaultCd<Undirected> > UndirectedDefaultCdToggle;



template<class Engine>
class CdSampler : public MetropolisHastings<Engine>{
protected:

	typedef std::set< std::pair<int,int> > DyadSet;
	DyadSet dyads;
	boost::shared_ptr< BinaryNet<Engine> > origNet;

public:

	CdSampler(Model<Engine> mod) : MetropolisHastings<Engine>(mod){
		boost::shared_ptr< DyadToggle<Engine, DefaultCd<Engine> > > tog(new DyadToggle<Engine, DefaultCd<Engine> >(*mod.network()));
		this->dyadToggle = tog;
		boost::shared_ptr< VertexToggle<Engine, DefaultVertex<Engine> > > tog1(
				new VertexToggle<Engine, DefaultVertex<Engine> >(*mod.network()));
		this->vertToggle = tog1;
		this->probDyad=.8;
	}
	CdSampler(Model<Engine> mod,int sampSize) : MetropolisHastings<Engine>(mod){
		Rcpp::List l;
		l.push_back(sampSize);
		boost::shared_ptr< DyadToggle<Engine, DefaultCd<Engine> > > tog(new DyadToggle<Engine, DefaultCd<Engine> >(l));
		tog->vSetNetwork(mod.network());
		this->dyadToggle = tog;
		boost::shared_ptr< VertexToggle<Engine, DefaultVertex<Engine> > > tog1(
				new VertexToggle<Engine, DefaultVertex<Engine> >(*mod.network()));
		this->vertToggle = tog1;
		this->probDyad=.8;
	}

	virtual ShallowCopyable* vShallowCopyUnsafe() const{
		return new CdSampler(*this);
	}


	void initialize(){
		MetropolisHastings<Engine>::initialize();
		dyads.clear();
		origNet = this->model->network()->clone();
	}

	void clearReferenceNetwork(){
		origNet = boost::shared_ptr< BinaryNet<Engine> >();
	}

	void setModel(const Model<Engine>& mod){
		MetropolisHastings<Engine>::setModel(mod);
		origNet = boost::shared_ptr< BinaryNet<Engine> >();
	}


	/*!
	 * Run MCMC for steps steps, then revert the network to its orig state
	 */
	double run(int steps){
		int accepted = 0;
		std::vector<std::pair<int,std::pair<int,double> > > contToggles;
		std::vector<std::pair<int,std::pair<int,int> > > disToggles;
		std::vector<std::pair<int,int> > tieToggles;

		double pd = this->probDyad;
		if(!this->model->hasAnyRandomVariables())
			pd = 1.0;
		if(!this->model->hasRandomGraph())
			pd = 0.0;

		bool isDyadToggle = false;
		for(int i=0;i<steps;i++){
			isDyadToggle = Rf_runif(0.0,1.0) < pd;
			if(isDyadToggle){
				this->dyadToggle->vGenerate();

				tieToggles = this->dyadToggle->vDyadToggles();
				contToggles.clear();
				disToggles.clear();
				//cout <<"toggles: " << tieToggles[0].first<<" "<<tieToggles[0].second<<"\n";
				//cout << "is tie: " << model->network()->hasEdge(tieToggles[0].first,tieToggles[0].second)<<"\n";
				double lastLik = this->model->logLik();
				for(int j=0;j<tieToggles.size();j++){
					this->model->dyadUpdate(tieToggles[j].first, tieToggles[j].second);
					//if(i<(tieToggles.size()-1))
					this->model->network()->toggle(tieToggles[j].first, tieToggles[j].second);
				}

				double lr = this->model->logLik() - lastLik;
				//cout << "lr:"<<lr<<"\n";
				lr += this->dyadToggle->vLogRatio();
				//cout << "v toggle lr:"<<dyadToggle.logRatio()<<"\n";
				if(lr>log(Rf_runif(0.0,1.0))){
					accepted++;
					this->dyadToggle->vTogglesAccepted(true);
					for(int j=0;j<tieToggles.size();j++)
						dyads.insert(tieToggles.at(j));
				}else{
					for(int j=0;j<tieToggles.size();j++){
						this->model->dyadUpdate(tieToggles[j].first, tieToggles[j].second);
						this->model->network()->toggle(tieToggles[j].first,tieToggles[j].second);
					}
					this->dyadToggle->vTogglesAccepted(false);
				}
			}else{
				Rf_error("no CD vertex toggling yet");
				this->vertToggle->vGenerate();
				tieToggles.clear();
				contToggles = this->vertToggle->vContVarChanges();
				disToggles = this->vertToggle->vDisVarChanges();
				if(contToggles.size()>0)
					::Rf_error("continuous variables unimplemented");
				double lastLik = this->model->logLik();

				std::vector<int> oldDisValues = std::vector<int>(disToggles.size(),-1);
				for(int j=0;j<disToggles.size();j++){
					oldDisValues[j] =this->model->network()->discreteVariableValue(
											disToggles[j].second.first,
											disToggles[j].first);
					//cout << "setting:"<<disToggles[j].first<<" "<<
					//		disToggles[j].second.first<<
					//		"\nwas:"<<oldDisValues.at(j)<<" to "<<
					//		disToggles.at(j).second.second<<"\n";
					this->model->discreteVertexUpdate(
							disToggles[j].first,
							disToggles[j].second.first ,
							disToggles[j].second.second);
					this->model->network()->setDiscreteVariableValue(
							disToggles[j].second.first,
							disToggles[j].first,
							disToggles[j].second.second);
				}


				double lr = this->model->logLik() - lastLik;
				//cout << "lr:"<<lr<<"\n";
				lr += this->vertToggle->vLogRatio();
				//cout << "d toggle lr:"<<vertToggle.logRatio()<<"\n";
				//cout <<" lr: "<<lr<<"\n";
				if(lr>log(Rf_runif(0.0,1.0))){
					accepted++;
					this->vertToggle->vTogglesAccepted(true);
				}else{
					//cout << "here"<<((int)disToggles.size())-1;
					for(int j=((int)disToggles.size())-1;j>=0;j--){
						//cout << "undoing change:"<<disToggles[j].first<<" "<<
						//		disToggles[j].second.first<<
						//		"\nto:"<<oldDisValues[j]<<"\n";
						this->model->discreteVertexUpdate(
							disToggles[j].first,
							disToggles[j].second.first ,
							oldDisValues.at(j));
						this->model->network()->setDiscreteVariableValue(
							disToggles[j].second.first,
							disToggles[j].first,
							oldDisValues.at(j));
					}
					this->vertToggle->vTogglesAccepted(false);
				}

			}
		}
		return ((double)accepted)/((double)steps);
	}

	void rollBackChanges(){
		for(DyadSet::iterator itr = dyads.begin(); itr != dyads.end(); ++itr){
			int from = (*itr).first;
			int to = (*itr).second;
			if(this->model->network()->hasEdge(from,to) != origNet->hasEdge(from,to)){
				this->model->dyadUpdate(from, to);
				this->model->network()->toggle(from,to);
			}
		}
		dyads.clear();
		MetropolisHastings<Engine>::initialize();
	}

	/*!
	 * Generates a list of BinaryNet objects using MCMC.
	 *
	 * Exposed to R
	 */
	Rcpp::List generateSample(int burnIn,int interval,int sampleSize){
		this->model->calculate();
		GetRNGstate();
		initialize();
		this->run(burnIn);
		Rcpp::List lis;
		for(int i=0;i<sampleSize-1;i++){
			R_CheckUserInterrupt();
			lis.push_back(this->model->network()->cloneR());
			rollBackChanges();
			this->run(interval);
		}
		lis.push_back(this->model->network()->cloneR());
		PutRNGstate();
		return lis;
	}

	/*!
	 * Generates a matrix of model statistics. Offset values are output as
	 * an attribute.
	 *
	 * Exposed to R
	 */
	NumericMatrix generateSampleStatistics(int burnIn,int interval,int sampleSize){
		std::vector<double> offs;
		std::vector<double> stats;
		this->model->calculate();
		NumericMatrix m(sampleSize,this->model->statistics().size());
		NumericMatrix off(sampleSize,this->model->offset().size());
		GetRNGstate();
		initialize();
		this->run(burnIn);
		for(int i=0;i<sampleSize;i++){
			R_CheckUserInterrupt();
			this->run(interval);
			stats = this->model->statistics();
			for(int j=0;j<stats.size();j++)
				m(i,j) = stats[j];
			offs = this->model->offset();
			for(int j=0;j<offs.size();j++)
				off(i,j) = offs[j];
			rollBackChanges();
		}
		PutRNGstate();
	    List lis;
	    lis.push_back(R_NilValue);
		lis.push_back(this->model->names());
		m.attr("dimnames") = lis;
		if(offs.size()>0)
			m.attr("offset") = off;
		return m;
	}

	/*!
	 * Generates a matrix of sample statistics. Additionally it applies
	 * the R Function supplimentalFunction to each sample and returns the
	 * result.
	 *
	 * Exported to R. Really shows the power of Rcpp
	 */
	List generateSampleStatisticsSupplimental(int burnIn, int interval,
			int sampleSize, Function supplimentalFunction){
		this->model->calculate();
		NumericMatrix m(sampleSize,this->model->statistics().size());
		List sup;
		GetRNGstate();
		initialize();
		this->run(burnIn);
		for(int i=0;i<sampleSize;i++){
			R_CheckUserInterrupt();
			this->run(interval);
			std::vector<double> stats = this->model->statistics();
			for(int j=0;j<stats.size();j++)
				m(i,j) = stats[j];
			sup.push_back(supplimentalFunction(*this->model->network()));
			rollBackChanges();
		}
		PutRNGstate();
		List lis;
		lis.push_back(m);
		lis.push_back(sup);
		return lis;
	}

	/** I guess I have to do this for the Rcpp Module...
	 *  copied from MetropolisHastings
	 */
	SEXP getModelR(){
		return wrap(*this->model);
	}
	void setDyadToggleType(std::string name, Rcpp::List params){
		MetropolisHastings<Engine>::setDyadToggleType(name,params);
		//this->dyadToggle = DyadTogglePtr(ToggleController<Engine>::getDyadToggle(name,params));
	}
	void setVertexToggleType(std::string name, Rcpp::List params){
		MetropolisHastings<Engine>::setVertexToggleType(name,params);
		//this->vertToggle = VertexTogglePtr(ToggleController<Engine>::getVertexToggle(name,params));
	}
	void setDyadProbability(double d){
		assert(d>=0.0 && d<=1.0);
		this->probDyad=d;
	}


};




/*!
 * A temporary experimental class to do Gibbsian CD
 */
template<class Engine>
class GibbsCdSampler : public MetropolisHastings<Engine>{
protected:

	typedef std::set< std::pair<int,int> > DyadSet;
	DyadSet dyads;
	boost::shared_ptr< BinaryNet<Engine> > origNet;

public:

	GibbsCdSampler(Model<Engine> mod) : MetropolisHastings<Engine>(mod){
		this->probDyad=.8;
	}
	GibbsCdSampler(Model<Engine> mod,int p) : MetropolisHastings<Engine>(mod,p){

	}

	virtual ShallowCopyable* vShallowCopyUnsafe() const{
		return new GibbsCdSampler(*this);
	}

	void initialize(){
		MetropolisHastings<Engine>::initialize();
		dyads.clear();
		origNet = this->model->network()->clone();
	}

	void clearReferenceNetwork(){
		origNet = boost::shared_ptr< BinaryNet<Engine> >();
	}

	void setModel(const Model<Engine>& mod){
		MetropolisHastings<Engine>::setModel(mod);
		origNet = boost::shared_ptr< BinaryNet<Engine> >();
	}


	/*!
	 * Run MCMC for steps steps, then revert the network to its orig state
	 */
	double run(int steps){
		int accepted = 0;
		std::vector<std::pair<int,std::pair<int,double> > > contToggles;
		std::vector<std::pair<int,std::pair<int,int> > > disToggles;
		std::vector<std::pair<int,int> > tieToggles;

		double pd = this->probDyad;
		if(!this->model->hasAnyRandomVariables())
			pd = 1.0;
		if(!this->model->hasRandomGraph())
			pd = 0.0;

		bool isDyadToggle = false;
		for(int i=0;i<steps;i++){
			isDyadToggle = Rf_runif(0.0,1.0) < pd;
			if(isDyadToggle){
				tieToggles = std::vector<std::pair<int,int> >(1,std::pair<int,int>(-1,-1));
				this->model->network()->randomDyad(tieToggles[0]);
				contToggles.clear();
				disToggles.clear();
				//cout <<"toggles: " << tieToggles[0].first<<" "<<tieToggles[0].second<<"\n";
				//cout << "is tie: " << model->network()->hasEdge(tieToggles[0].first,tieToggles[0].second)<<"\n";
				double lastLik = this->model->logLik();
				for(int j=0;j<tieToggles.size();j++){
					this->model->dyadUpdate(tieToggles[j].first, tieToggles[j].second);
					//if(i<(tieToggles.size()-1))
					this->model->network()->toggle(tieToggles[j].first, tieToggles[j].second);
				}

				//Gibbs prob of change
				double p = 1.0 / ( 1.0 + exp(lastLik - this->model->logLik()));

				if( Rf_runif(0.0,1.0) <= p ){
					accepted++;
					for(int j=0;j<tieToggles.size();j++)
						dyads.insert(tieToggles.at(j));
				}else{
					for(int j=0;j<tieToggles.size();j++){
						this->model->dyadUpdate(tieToggles[j].first, tieToggles[j].second);
						this->model->network()->toggle(tieToggles[j].first,tieToggles[j].second);
					}
				}
			}else{
				Rf_error("no CD vertex toggling yet");
			}
		}
		return ((double)accepted)/((double)steps);
	}

	void rollBackChanges(){
		for(DyadSet::iterator itr = dyads.begin(); itr != dyads.end(); ++itr){
			int from = (*itr).first;
			int to = (*itr).second;
			if(this->model->network()->hasEdge(from,to) != origNet->hasEdge(from,to)){
				this->model->dyadUpdate(from, to);
				this->model->network()->toggle(from,to);
			}
		}
		dyads.clear();
		MetropolisHastings<Engine>::initialize();
	}

	/*!
	 * Generates a list of BinaryNet objects using MCMC.
	 *
	 * Exposed to R
	 */
	Rcpp::List generateSample(int burnIn,int interval,int sampleSize){
		this->model->calculate();
		GetRNGstate();
		initialize();
		this->run(burnIn);
		Rcpp::List lis;
		for(int i=0;i<sampleSize-1;i++){
			R_CheckUserInterrupt();
			lis.push_back(this->model->network()->cloneR());
			rollBackChanges();
			this->run(interval);
		}
		lis.push_back(this->model->network()->cloneR());
		PutRNGstate();
		return lis;
	}

	/*!
	 * Generates a matrix of model statistics. Offset values are output as
	 * an attribute.
	 *
	 * Exposed to R
	 */
	NumericMatrix generateSampleStatistics(int burnIn,int interval,int sampleSize){
		std::vector<double> offs;
		std::vector<double> stats;
		this->model->calculate();
		NumericMatrix m(sampleSize,this->model->statistics().size());
		NumericMatrix off(sampleSize,this->model->offset().size());
		GetRNGstate();
		initialize();
		this->run(burnIn);
		for(int i=0;i<sampleSize;i++){
			R_CheckUserInterrupt();
			this->run(interval);
			stats = this->model->statistics();
			for(int j=0;j<stats.size();j++)
				m(i,j) = stats[j];
			offs = this->model->offset();
			for(int j=0;j<offs.size();j++)
				off(i,j) = offs[j];
			rollBackChanges();
		}
		PutRNGstate();
	    List lis;
	    lis.push_back(R_NilValue);
		lis.push_back(this->model->names());
		m.attr("dimnames") = lis;
		if(offs.size()>0)
			m.attr("offset") = off;
		return m;
	}

	/*!
	 * Generates a matrix of sample statistics. Additionally it applies
	 * the R Function supplimentalFunction to each sample and returns the
	 * result.
	 *
	 * Exported to R. Really shows the power of Rcpp
	 */
	List generateSampleStatisticsSupplimental(int burnIn, int interval,
			int sampleSize, Function supplimentalFunction){
		this->model->calculate();
		NumericMatrix m(sampleSize,this->model->statistics().size());
		List sup;
		GetRNGstate();
		initialize();
		this->run(burnIn);
		for(int i=0;i<sampleSize;i++){
			R_CheckUserInterrupt();
			this->run(interval);
			std::vector<double> stats = this->model->statistics();
			for(int j=0;j<stats.size();j++)
				m(i,j) = stats[j];
			sup.push_back(supplimentalFunction(*this->model->network()));
			rollBackChanges();
		}
		PutRNGstate();
		List lis;
		lis.push_back(m);
		lis.push_back(sup);
		return lis;
	}

	/** I guess I have to do this for the Rcpp Module...
	 *  copied from MetropolisHastings
	 */
	SEXP getModelR(){
		return wrap(*this->model);
	}
	void setDyadToggleType(std::string name, Rcpp::List params){
		MetropolisHastings<Engine>::setDyadToggleType(name,params);
		//this->dyadToggle = DyadTogglePtr(ToggleController<Engine>::getDyadToggle(name,params));
	}
	void setVertexToggleType(std::string name, Rcpp::List params){
		MetropolisHastings<Engine>::setVertexToggleType(name,params);
		//this->vertToggle = VertexTogglePtr(ToggleController<Engine>::getVertexToggle(name,params));
	}
	void setDyadProbability(double d){
		assert(d>=0.0 && d<=1.0);
		this->probDyad=d;
	}


};




/*!
 * A temporary experimental class to do Gibbsian CD
 */
template<class Engine>
class GibbsCdSampler2 : public MetropolisHastings<Engine>{
protected:

	typedef std::set< std::pair<int,int> > DyadSet;
	DyadSet dyads;
	boost::shared_ptr< BinaryNet<Engine> > origNet;
	int ss;
public:
	GibbsCdSampler2(Model<Engine> mod) : MetropolisHastings<Engine>(mod){
		this->probDyad=.8;
		ss = mod.network()->size();
	}
	GibbsCdSampler2(Model<Engine> mod,int s) : MetropolisHastings<Engine>(mod){
		this->probDyad=.8;
		ss = s;
	}
	GibbsCdSampler2(Model<Engine> mod,double p) : MetropolisHastings<Engine>(mod,p){
		ss = mod.network()->size();
	}

	virtual ShallowCopyable* vShallowCopyUnsafe() const{
		return new GibbsCdSampler2(*this);
	}

	void initialize(){
		MetropolisHastings<Engine>::initialize();
		dyads.clear();
		origNet = this->model->network()->clone();
	}

	void clearReferenceNetwork(){
		origNet = boost::shared_ptr< BinaryNet<Engine> >();
	}

	void setModel(const Model<Engine>& mod){
		MetropolisHastings<Engine>::setModel(mod);
		origNet = boost::shared_ptr< BinaryNet<Engine> >();
	}


	/*!
	 * Run MCMC for steps steps, then revert the network to its orig state
	 */
	double run(int steps){
		std::vector<int> order;
		sampleWithoutReplacement(this->model->network()->size(),ss,order);
		int n1 = floor(Rf_runif(0.0,ss));
		int node = order.at(n1);
		int accepted = 0;
		std::vector<std::pair<int,std::pair<int,double> > > contToggles;
		std::vector<std::pair<int,std::pair<int,int> > > disToggles;
		std::vector<std::pair<int,int> > tieToggles;

		double pd = this->probDyad;
		if(!this->model->hasAnyRandomVariables())
			pd = 1.0;
		if(!this->model->hasRandomGraph())
			pd = 0.0;

		bool isDyadToggle = false;
		for(int i=0;i<steps;i++){
			isDyadToggle = Rf_runif(0.0,1.0) < pd;
			if(isDyadToggle){
				int nbr = floor(Rf_runif(0.0,ss-1.0));
				if(nbr>=n1)
					nbr++;
				nbr = order.at(nbr);
				tieToggles = std::vector<std::pair<int,int> >(1,std::pair<int,int>(node,nbr));
				//this->model->network()->randomDyad(tieToggles[0]);
				contToggles.clear();
				disToggles.clear();
				//cout <<"toggles: " << tieToggles[0].first<<" "<<tieToggles[0].second<<"\n";
				//cout << "is tie: " << model->network()->hasEdge(tieToggles[0].first,tieToggles[0].second)<<"\n";
				double lastLik = this->model->logLik();
				for(int j=0;j<tieToggles.size();j++){
					this->model->dyadUpdate(tieToggles[j].first, tieToggles[j].second);
					//if(i<(tieToggles.size()-1))
					this->model->network()->toggle(tieToggles[j].first, tieToggles[j].second);
				}

				//Gibbs prob of change
				double p = 1.0 / ( 1.0 + exp(lastLik - this->model->logLik()));

				if( Rf_runif(0.0,1.0) <= p ){
					accepted++;
					for(int j=0;j<tieToggles.size();j++)
						dyads.insert(tieToggles.at(j));
				}else{
					for(int j=0;j<tieToggles.size();j++){
						this->model->dyadUpdate(tieToggles[j].first, tieToggles[j].second);
						this->model->network()->toggle(tieToggles[j].first,tieToggles[j].second);
					}
				}
			}else{
				Rf_error("no CD vertex toggling yet");
			}
		}
		return ((double)accepted)/((double)steps);
	}

	void rollBackChanges(){
		for(DyadSet::iterator itr = dyads.begin(); itr != dyads.end(); ++itr){
			int from = (*itr).first;
			int to = (*itr).second;
			if(this->model->network()->hasEdge(from,to) != origNet->hasEdge(from,to)){
				this->model->dyadUpdate(from, to);
				this->model->network()->toggle(from,to);
			}
		}
		dyads.clear();
		MetropolisHastings<Engine>::initialize();
	}

	/*!
	 * Generates a list of BinaryNet objects using MCMC.
	 *
	 * Exposed to R
	 */
	Rcpp::List generateSample(int burnIn,int interval,int sampleSize){
		this->model->calculate();
		GetRNGstate();
		initialize();
		this->run(burnIn);
		Rcpp::List lis;
		for(int i=0;i<sampleSize-1;i++){
			R_CheckUserInterrupt();
			lis.push_back(this->model->network()->cloneR());
			rollBackChanges();
			this->run(interval);
		}
		lis.push_back(this->model->network()->cloneR());
		PutRNGstate();
		return lis;
	}

	/*!
	 * Generates a matrix of model statistics. Offset values are output as
	 * an attribute.
	 *
	 * Exposed to R
	 */
	NumericMatrix generateSampleStatistics(int burnIn,int interval,int sampleSize){
		std::vector<double> offs;
		std::vector<double> stats;
		this->model->calculate();
		NumericMatrix m(sampleSize,this->model->statistics().size());
		NumericMatrix off(sampleSize,this->model->offset().size());
		GetRNGstate();
		initialize();
		this->run(burnIn);
		for(int i=0;i<sampleSize;i++){
			R_CheckUserInterrupt();
			this->run(interval);
			stats = this->model->statistics();
			for(int j=0;j<stats.size();j++)
				m(i,j) = stats[j];
			offs = this->model->offset();
			for(int j=0;j<offs.size();j++)
				off(i,j) = offs[j];
			rollBackChanges();
		}
		PutRNGstate();
	    List lis;
	    lis.push_back(R_NilValue);
		lis.push_back(this->model->names());
		m.attr("dimnames") = lis;
		if(offs.size()>0)
			m.attr("offset") = off;
		return m;
	}

	/*!
	 * Generates a matrix of sample statistics. Additionally it applies
	 * the R Function supplimentalFunction to each sample and returns the
	 * result.
	 *
	 * Exported to R. Really shows the power of Rcpp
	 */
	List generateSampleStatisticsSupplimental(int burnIn, int interval,
			int sampleSize, Function supplimentalFunction){
		this->model->calculate();
		NumericMatrix m(sampleSize,this->model->statistics().size());
		List sup;
		GetRNGstate();
		initialize();
		this->run(burnIn);
		for(int i=0;i<sampleSize;i++){
			R_CheckUserInterrupt();
			this->run(interval);
			std::vector<double> stats = this->model->statistics();
			for(int j=0;j<stats.size();j++)
				m(i,j) = stats[j];
			sup.push_back(supplimentalFunction(*this->model->network()));
			rollBackChanges();
		}
		PutRNGstate();
		List lis;
		lis.push_back(m);
		lis.push_back(sup);
		return lis;
	}

	/** I guess I have to do this for the Rcpp Module...
	 *  copied from MetropolisHastings
	 */
	SEXP getModelR(){
		return wrap(*this->model);
	}
	void setDyadToggleType(std::string name, Rcpp::List params){
		MetropolisHastings<Engine>::setDyadToggleType(name,params);
		//this->dyadToggle = DyadTogglePtr(ToggleController<Engine>::getDyadToggle(name,params));
	}
	void setVertexToggleType(std::string name, Rcpp::List params){
		MetropolisHastings<Engine>::setVertexToggleType(name,params);
		//this->vertToggle = VertexTogglePtr(ToggleController<Engine>::getVertexToggle(name,params));
	}
	void setDyadProbability(double d){
		assert(d>=0.0 && d<=1.0);
		this->probDyad=d;
	}


};



} /* namespace ernm */
#endif /* CDSAMPLER_H_ */
