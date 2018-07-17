/**
 ** This software is covered by the "BSD 2-Clause License"
 ** and any user of this software or source file is bound by the terms therein.
 ** 
 ** Redistribution and use in source and binary forms, with or without
 ** modification, are permitted provided that the following conditions are met:
 **
 ** - Redistributions of source code must retain the above copyright notice, this
 **   list of conditions and the following disclaimer.
 **
 ** - Redistributions in binary form must reproduce the above copyright notice,
 **   this list of conditions and the following disclaimer in the documentation
 **   and/or other materials provided with the distribution.
 **
 **
 ** This software is an implementation of the algorithm described in:
 ** The colored longest common prefix array computed via sequential scans
 ** SPIRE 2018
 ** by F. Garofalo, G. Rosone, M. Sciortino and D. Verzotto
 ** 
 ** 
 ** Supported by the project Italian MIUR-SIR CMACBioSeq 
 ** (``Combinatorial methods for analysis and compression of biological sequences'') 
 ** grant n.~RBSI146R5L.
 ** 
 ** 
 ** Copyright by the above authors.
 ** 
 **
 ** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 ** AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 ** IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 ** DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 ** FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 ** DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 ** SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 ** CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 ** OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 ** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 **/
#include "MultiACS.h"
#include "Parameters.h"
#include "Types.h"
#include "Tools.h"
#include "Writer.h"
#include "Reader.h"
#include "StackedDGenerator.h"
#include "GESAConverter.h"
#include "malloc_count/malloc_count.h"
#include <vector>
#include <unistd.h>
#include <sstream>
#include <string>
#include <cstring>
#include <cmath>



using namespace std;

namespace multi_acs {

MultiACS::MultiACS(MultiACSParameters* params) :
	params(params),
	global_collection(params->target_collection_file_name, true, params->file_format, params->verbose),
	reference_color(params->reference_color) {

	if(global_collection.colors.count(params->reference_color) == 0) {
		ostringstream err_message;
		err_message << "Couldn't find reference color in target collection";
		Error::stopWithError(C_MultiACS_ClassName, __func__, err_message.str());
	}
	global_collection.printCollectionInfo();

	generateD();
	generateLCP();

}

MultiACS::MultiACS(MultiACSParameters* params, CollectionInfo &collection) :
	params(params),
	global_collection(collection),
	reference_color(params->reference_color) { }

MultiACS::~MultiACS() {
	// TODO Auto-generated destructor stub
}

void MultiACS::generateD() {
	StackedDGenerator d_gen(params);
	d_gen.generateD();
}

void MultiACS::generateLCP() {
	// TODO NOT YET IMPLEMENTED
}

void MultiACS::computeACS() {

	SequenceNumber m = global_collection.getSequenceNumber();
	SequenceLength n_x = global_collection.getSequenceLength(reference_color);

	LetterNumber score_x[m]{0};
	memset (score_x, 0, m);
	LetterNumber score_r[m]{0};
	memset (score_r, 0, m);

	clock_t start = clock();
	time_t start_wc = time(NULL);
	double elapsed, elapsed_wc;

	cout << "ACS Computation Starting\n";
	cout << "Current memory: " << malloc_count_current() << " bytes" << endl;

	forwardComputation(score_r);
	backwardComputation(score_x);

	//cout << "ACS distance_file_name\n";
	
	FileName distance_file_name(params->output_file_name, C_DistanceFileExt);
	FILE* distance_file = fopen(distance_file_name.c_str(), "w");
	if(distance_file == nullptr) {
		ostringstream err_message;
		err_message << "Couldn't open file " << distance_file_name.str();
		Error::stopWithError(C_MultiACS_ClassName, __func__, err_message.str());
	}

//	AlphabetSize sigma = global_collection.getAlphabetSize() - 1;
//	double base = log(sigma);
	SequenceLength s1, s2;
	DistanceValue sumS1, sumS2;
	DistanceValue distance = 0;
	for(SequenceNumber r = 0; r < m; ++r) {
//		cout << "SequenceNumber: " << r << "\n";
		if(r != reference_color) {
			SequenceLength n_r = global_collection.getSequenceLength(r);
//			distance = (n_x - 1)*(log(n_r - 1)/base)/(2*score_x[r])
//					- (log(n_x - 1)/base)/(n_x)
//					+ (n_r - 1)*(log(n_x - 1)/base)/(2*score_r[r])
//					- (log(n_r - 1)/base)/(n_r);
			s1 = n_x - 1;
			s2 = n_r - 1;
			sumS1 = score_x[r];
			sumS2 = score_r[r];

			if(params->verbose) {
				cout << "|Seq" << reference_color << "|=" << n_x - 1 <<
						"\t|Seq" << r << "|=" << n_r - 1 << '\n';
				cout << "score(" << reference_color << "," << r << ")="
						<< score_x[r] <<
						"\tscore(" << r << "," << reference_color << ")="
						<< score_r[r] << endl;
			}

			distance = ((log10(s1)/(sumS2/s2))-((2.0*log10(s2))/s2) + (log10(s2)/(sumS1/s1))-((2.0*log10(s1))/s1))*0.5;

			fprintf(distance_file, "%f\t", distance);
		}
		else {
			fprintf(distance_file, "0\t");
		}
	}

	fclose(distance_file);

	cout << "ACS Computation End\n";
	cout << "Current memory: " << malloc_count_current() << " bytes" << endl;

	elapsed = (clock()-start)/(double)(CLOCKS_PER_SEC);
	elapsed_wc = difftime(time(NULL),start_wc);

	cout << "### ACS Computation Reporting ###\n";
	cout << "Elapsed time: " << elapsed << " secs\n";
	cout << "Wall Clock time: " << elapsed_wc << " secs\n";
	cout << "Peak memory: " << malloc_count_peak() << " bytes\n";
	cout << endl;

}

void MultiACS::forwardComputation(LetterNumber score_r[]) {

	const SequenceNumber m = global_collection.getSequenceNumber();
	const SequenceLength n_x = global_collection.getSequenceLength(reference_color);

	FileName id_file_name(params->target_collection_file_name, C_IdFileExt);
	FILE* id_file = fopen(id_file_name.c_str(), "rb");
	if(id_file == nullptr) {
		ostringstream err_message;
		err_message << "Couldn't open file " << id_file_name.str();
		Error::stopWithError(C_MultiACS_ClassName, __func__, err_message.str());
	}
//	cout << "Reading from" << id_file_name.str() << endl;

	FileName lcp_file_name(params->target_collection_file_name, C_LcpFileExt);
	FILE* lcp_file = fopen(lcp_file_name.c_str(), "rb");
	if(lcp_file == nullptr) {
		ostringstream err_message;
		err_message << "Couldn't open file " << lcp_file_name.str();
		Error::stopWithError(C_MultiACS_ClassName, __func__, err_message.str());
	}
//	cout << "Reading from" << lcp_file_name.str() << endl;

	FileName d_file_name(params->output_file_name, C_DynBlockFileExt);
	FILE* d_file = fopen(d_file_name.c_str(), "rb");
	if(d_file == nullptr) {
		ostringstream err_message;
		err_message << "Couldn't open file " << d_file_name.str();
		Error::stopWithError(C_MultiACS_ClassName, __func__, err_message.str());
	}
//	cout << "Reading from" << d_file_name.str() << endl;

	FileName lcp_x_file_name = FileName(params->reference_sequence_file_name, C_LcpFileExt);
	FILE* lcp_x_file = fopen(lcp_x_file_name.c_str(), "rb");
	if(lcp_x_file == nullptr) {
		ostringstream err_message;
		err_message << "Couldn't open file " << lcp_x_file_name.str();
		Error::stopWithError(C_MultiACS_ClassName, __func__, err_message.str());
	}
//	cout << "Reading from" << lcp_x_file_name.str() << endl;

	FileName cLCP_x_file_name(params->output_file_name, C_PartialCLcpFileExt);
	FILE* cLCP_x_file = fopen(cLCP_x_file_name.c_str(), "wb");
	if(cLCP_x_file == nullptr) {
		ostringstream err_message;
		err_message << "Couldn't open file " << cLCP_x_file_name.str();
		Error::stopWithError(C_MultiACS_ClassName, __func__, err_message.str());
	}
//	cout << "Writing from" << cLCP_x_file_name.str() << endl;

	AllocableMemory A = params->memory_amount/sizeof(SequenceLength);

	SequenceLength Q = ceil((double) A/m);
	Q = Q > n_x ? n_x : Q;

	SequenceLength cLCP_x[Q + 1][m];
	memset(cLCP_x, 0, sizeof(cLCP_x));
	bool LcLCPbit_x[m];

	IdReader id_reader(id_file);
	LCPReader lcp_reader(lcp_file);
	LCPReader lcp_x_reader(lcp_x_file);
	DReader d_reader(d_file);

	SequenceLength h_x = 0;
	SequenceLength h_x_idx = 0;
	SequenceLength alpha = C_MaxSequenceLength;
	SequenceLength k = 0;
	SequenceNumber id;
	SequenceLength lcp_x_value;

	cout << "cLCP Forward Computation\n";
	cout << "Current memory: " << malloc_count_current() << " bytes" << endl;

	lcp_x_reader.readLCPValue(lcp_x_value);
	while(id_reader.readSequenceId(id)) {

		SequenceLength lcp_value, d_value;
		lcp_reader.readLCPValue(lcp_value);
		d_reader.readDValue(d_value);

		if(id != reference_color) {

			if(lcp_value < alpha)
				alpha = lcp_value;

			if((d_value > 0) && (d_value - 1 > k))
				k = d_value - 1;

			if(h_x == 0) {
				score_r[id] += k;
				cLCP_x[h_x_idx][id] = k;
				if(params->verbose) {
					cout << id << ": " << k << endl;
				}
//				cout << id << ": " << k << endl;
			}
			else if(alpha > lcp_x_value) {
				score_r[id] += alpha;
				if(h_x < n_x) {
					cLCP_x[h_x_idx][id] = lcp_x_value;
				}
				if(params->verbose) {
					cout << id << ": " << alpha << endl;
				}
//				cout << id << ": " << alpha << endl;
			}
			else {
				score_r[id] += max(max(alpha,k),lcp_x_value);
				cLCP_x[h_x_idx][id] = max(k, lcp_x_value);
				if(params->verbose) {
					cout << id << ": " << max(max(alpha,k),lcp_x_value) << endl;
				}
//				cout << id << ": " << max(max(alpha,k),lcp_x_value) << endl;
			}

			if(h_x > 0 && LcLCPbit_x[id] == false) {
				if(alpha > cLCP_x[h_x_idx - 1][id])
					cLCP_x[h_x_idx - 1][id] = alpha;
				LcLCPbit_x[id] = true;
			}
		}
		else {
			if(h_x > 0) {
				for(SequenceNumber r = 0; r < m; ++r) {
					cLCP_x[h_x_idx][r] = max(min(cLCP_x[h_x_idx - 1][r], lcp_x_value), cLCP_x[h_x_idx][r]);
				}

				if(params->verbose) {
					cout << "[" << h_x - 1 << "]:";
					for(SequenceNumber r = 0; r < m; ++r) {
						cout << " " << cLCP_x[h_x_idx - 1][r];
					}
					cout << endl;
				}
//				cout << "[" << h_x - 1 << "]:";
//				for(SequenceNumber r = 0; r < m; ++r) {
//					cout << " " << cLCP_x[h_x_idx - 1][r];
//				}
//				cout << endl;
			}

			++h_x;
			++h_x_idx;
			if(h_x_idx == Q + 1) {
				fwrite(cLCP_x[0], sizeof(SequenceLength), Q*m, cLCP_x_file);
				for(SequenceNumber r = 0; r < m; ++r) {
					cLCP_x[0][r] = cLCP_x[Q][r];
					for(SequenceLength j = 1; j < Q + 1; ++j){
						cLCP_x[j][r] = 0;
					}
				}
				h_x_idx = 1;
			}
			alpha = C_MaxSequenceLength;
			k = 0;
			for(SequenceNumber r = 0; r < m; ++r) {
				LcLCPbit_x[r] = false;
			}

			if(!lcp_x_reader.readLCPValue(lcp_x_value))
				lcp_x_value = 0;
		}
	}

	if(params->verbose) {
		cout << "[" << h_x - 1 << "]:";
		for(SequenceNumber r = 0; r < m; ++r) {
			cout << " " << cLCP_x[h_x_idx - 1][r];
		}
		cout << endl;
	}
//	cout << "[" << h_x - 1 << "]:";
//	for(SequenceNumber r = 0; r < m; ++r) {
//		cout << " " << cLCP_x[h_x_idx - 1][r];
//	}
//	cout << endl;

	fwrite(cLCP_x[0], sizeof(SequenceLength), h_x_idx*m, cLCP_x_file);

	fclose(cLCP_x_file);
	fclose(lcp_x_file);
	fclose(d_file);
	fclose(lcp_file);
	fclose(id_file);
}

void MultiACS::backwardComputation(LetterNumber score_x[]) {

	SequenceNumber m = global_collection.getSequenceNumber();
	SequenceLength n_x = global_collection.getSequenceLength(reference_color);

	FileName lcp_x_file_name(params->reference_sequence_file_name, C_LcpFileExt);
	FILE* lcp_x_file = fopen(lcp_x_file_name.c_str(), "rb");
	if(lcp_x_file == nullptr) {
		ostringstream err_message;
		err_message << "Couldn't open file " << lcp_x_file_name.str();
		Error::stopWithError(C_MultiACS_ClassName, __func__, err_message.str());
	}
//	cout << "Reading from" << lcp_x_file_name.c_str() << endl;

	FileName cLCP_x_file_name(params->output_file_name, C_PartialCLcpFileExt);
	FILE* cLCP_x_file = fopen(cLCP_x_file_name.c_str(), "rb");
	if(cLCP_x_file == nullptr) {
		ostringstream err_message;
		err_message << "Couldn't open file " << cLCP_x_file_name.str();
		Error::stopWithError(C_MultiACS_ClassName, __func__, err_message.str());
	}
//	cout << "Reading from" << cLCP_x_file_name.str() << endl;

	AllocableMemory A = params->memory_amount/sizeof(SequenceLength);

	SequenceLength Q = ceil((double) A/m);
	Q = Q > n_x ? n_x : Q;

	//SequenceLength cLCP_x[Q + 1][m];
	vector<vector<SequenceLength>> cLCP_x(Q + 1, vector<SequenceLength>(m, 0));
	SequenceLength lcp_x[Q + 1];
	SequenceLength H = ceil((double) n_x/Q);
	SequenceLength q, e;
	q = n_x;

	cout << "cLCP Backward Computation\n";
	cout << "Current memory: " << malloc_count_current() << " bytes" << endl;


	for(SequenceLength h = 1; h <= H ; ++h) {
		//if (h % 1000 == 0)
			//cout << "h = " << h << "\n";
		if(Q > q)
			Q = q;
		q = q - Q;
		// e = (h != H) ? 1 : 0;
		e = (h != 1) ? 1 : 0;
		fseek(cLCP_x_file, q*m*sizeof(SequenceLength), SEEK_SET);
		fseek(lcp_x_file, q*sizeof(SequenceLength), SEEK_SET);
		for(SequenceLength j = 0; j < Q + e; ++j)
			fread(cLCP_x[j].data(), sizeof(SequenceLength), m, cLCP_x_file);
		fread(lcp_x, sizeof(SequenceLength), Q + e, lcp_x_file);
		for(SequenceLength k = Q; k > 0 ; --k) {
			for(SequenceNumber r = 0; r < m; ++r) {
				if(q + k == n_x) {
					score_x[r] += cLCP_x[k - 1][r];
				}
				else {
					cLCP_x[k - 1][r] = max(min(cLCP_x[k][r], lcp_x[k]), cLCP_x[k - 1][r]);
					score_x[r] += cLCP_x[k - 1][r];
				}
			}
			if(params->verbose) {
				cout << "[" << q + k - 1 << "]:";
				for(SequenceNumber r = 0; r < m; ++r) {
					cout << " " << cLCP_x[k - 1][r];
				}
				cout << endl;
			}
		}

	}

	//delete [] cLCP_x;
	cout << "END -------- cLCP Backward Computation\n";
	
	fclose(cLCP_x_file);
	fclose(lcp_x_file);

}



} /* namespace multi_acs */

void printUsage() {
	cout << "Usage: [-h] [-v] [-f input_format] [-Q amount] ref_seq target_seqs ref_color output" << endl;
}

using namespace multi_acs;

int main(int argc, char* argv[]) {

	bool verbose = false;
	int input_format = 0;
	string reference_seq_file_name, target_collection_file_name, output_file_name;
	SequenceNumber reference_color;
	AllocableMemory memory_amount = BUFFER_SIZE*sizeof(SequenceLength);

	int o;
	while((o = getopt(argc, argv, "vhf:Q:")) != -1) {
		switch(o) {
			case 'v':
				verbose = true;
				break;
			case 'f':
				input_format = atoi(optarg);
				break;
			case 'Q':
				memory_amount = atoi(optarg);
				break;
			case 'h':
			default:
				printUsage();
				exit(EXIT_FAILURE);
		}
	}

	if(optind == argc - 4) {
		reference_seq_file_name = string(argv[optind++]);
	}
	else {
		printUsage();
		ostringstream err_message;
		err_message << "Missing Reference Sequence file name";
		Error::stopWithError(C_MultiACS_ClassName, __func__, err_message.str());
	}

	if(optind == argc - 3) {
		target_collection_file_name = string(argv[optind++]);
	}
	else {
		printUsage();
		ostringstream err_message;
		err_message << "Missing target collection file name";
		Error::stopWithError(C_MultiACS_ClassName, __func__, err_message.str());
	}

	if(optind == argc - 2) {
		reference_color = atoi(argv[optind++]);
	}
	else {
		printUsage();
		ostringstream err_message;
		err_message << "Missing reference color";
		Error::stopWithError(C_MultiACS_ClassName, __func__, err_message.str());
	}

	if(optind == argc - 1) {
		output_file_name = string(argv[optind++]);
	}
	else {
		printUsage();
		ostringstream err_message;
		err_message << "Missing output file name";
		Error::stopWithError(C_MultiACS_ClassName, __func__, err_message.str());
	}

	MultiACSParameters* params = new MultiACSParameters(verbose,
			reference_seq_file_name,
			target_collection_file_name,
			input_format,
			reference_color,
			output_file_name, memory_amount);
	params->printParameters();

	/*  Build separated files from GESA file	
	if(input_format == 1) {
		GESAConverter::extractFromGESA(reference_seq_file_name);
		GESAConverter::extractFromGESA(target_collection_file_name);
	}
	*/
	
/*
	CollectionInfo reference_sequence(reference_seq_file_name,
		0,
		false,// do not collect symbols/colors frequencies
		params->verbose);
*/		


/* Find info from gesa   
	CollectionInfo collection(params->target_collection_file_name,
			0,
			true,
			params->verbose);
	if(collection.colors.count(params->reference_color) == 0) {
		ostringstream err_message;
		err_message << "Couldn't find reference color in target collection";
		Error::stopWithError(C_MultiACS_ClassName, __func__, err_message.str());
	}
	collection.printCollectionInfo();
*/

//Take info from file if no pre-processing 
//(in this case, the above lines are not necessary)
//	CollectionInfo collection(params->target_collection_file_name);
//	collection.printCollectionInfo();

//Load the length if there exists a file containing the lengths
CollectionInfo collection;
collection.loadCollectionLengths(params->target_collection_file_name);
collection.printCollectionInfo();

	clock_t start = clock();
	time_t start_wc = time(NULL);
	double elapsed, elapsed_wc;

	cout << "Global Computation Starting\n";
	cout << "Current memory: " << malloc_count_current() << " bytes" << endl;

/**/
	StackedDGenerator d_gen(params);
	d_gen.generateD();
/**/

	MultiACS acs(params, collection);
	acs.computeACS();

	cout << "Global Computation End\n";
	cout << "Current memory: " << malloc_count_current() << " bytes" << endl;

	elapsed = (clock()-start)/(double)(CLOCKS_PER_SEC);
	elapsed_wc = difftime(time(NULL),start_wc);

	cout << "### Global Computation Reporting ###\n";
	cout << "Elapsed time: " << elapsed << " secs\n";
	cout << "Wall Clock time: " << elapsed_wc << " secs\n";
	cout << "Peak memory: " << malloc_count_peak() << " bytes\n";
	cout << endl;

}
