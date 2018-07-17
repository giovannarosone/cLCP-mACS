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
#include "StackedDGenerator.h"
#include "Types.h"
#include "Tools.h"
#include "Reader.h"
#include "Writer.h"
#include "malloc_count/malloc_count.h"
#include <sstream>
#include <list>

using namespace std;

namespace multi_acs {

StackedDGenerator::StackedDGenerator(MultiACSParameters* params) :
	params(params) { }

StackedDGenerator::~StackedDGenerator() {
	// TODO Auto-generated destructor stub
}

bool StackedDGenerator::mapColor(SequenceNumber id) {
	return (id == params->reference_color);
}

void StackedDGenerator::generateDPairs() {

	std::string input_file_name = params->target_collection_file_name + C_GESAExt;
	FILE* f_ESA = fopen(input_file_name.c_str(), "rb");
	if(f_ESA == nullptr) {
		std::ostringstream err_message;
		err_message << "Couldn't open file " << input_file_name;
		Error::stopWithError(C_StackedDGenerator_ClassName, __func__, err_message.str());
	}
//	std::cout << "Read from EGSA File: " << input_file_name << std::endl;

	std::string d_file_name = params->output_file_name + C_DynBlockFileExt;
	FILE* d_file = fopen(d_file_name.c_str(), "wb");
	if(d_file == nullptr) {
		std::ostringstream err_message;
		err_message << "Couldn't open file " << d_file_name;
		Error::stopWithError(C_StackedDGenerator_ClassName, __func__, err_message.str());
	}
//	std::cout << "Write to D File: " << d_file_name << std::endl;

	list<stackedLCPInterval> stacked_list;
	stackedLCPInterval lcp_interval;
	t_GSA buffer[BUFFER_SIZE];
	LetterNumber num_read = 0, k = 0;
	SequenceLength top_lcp = 0, max_common_lcp = 0;
	bool current_color, successive_color;

	clock_t start = clock();
	time_t start_wc = time(NULL);
	double elapsed, elapsed_wc;

	cout << "D Computation Starting\n";
	cout << "Current memory: " << malloc_count_current() << " bytes" << endl;

	if((num_read = fread(buffer, sizeof(buffer[0]), BUFFER_SIZE, f_ESA)) > 0) {

		while(buffer[k].lcp == 0) ++k;

		lcp_interval = {
				k - 1,				// pos
				buffer[k].lcp,		// lcp value
		};
		stacked_list.push_front(lcp_interval);
		// cout << stacked_list.size() << endl;
		top_lcp = buffer[k].lcp;
		current_color = mapColor(buffer[k].text);
		++k;

		for(LetterNumber i = k; i < num_read; ++i, ++k) {

			successive_color = mapColor(buffer[i].text);

			if(buffer[i].lcp != 0) {

				if(buffer[i].lcp > top_lcp) {

					lcp_interval = {
							k - 1,				// pos
							buffer[i].lcp,		// lcp value
					};
					stacked_list.push_front(lcp_interval);
					// cout << stacked_list.size() << endl;
					top_lcp = buffer[i].lcp;

				} else if (buffer[i].lcp < top_lcp) {

					LetterNumber ini_pos = i - 1;
					while(!stacked_list.empty() && buffer[i].lcp < top_lcp) {
						ini_pos = stacked_list.front().pos;
						stacked_list.pop_front();
						// cout << stacked_list.size() << endl;
						if(!stacked_list.empty())
							top_lcp = stacked_list.front().lcp;
						else
							top_lcp = 0;
					}

					if(buffer[i].lcp > max_common_lcp) {

						if(buffer[i].lcp > top_lcp) {

							lcp_interval = {
									ini_pos,			// pos
									buffer[i].lcp,	// lcp value
							};
							stacked_list.push_front(lcp_interval);
							// cout << stacked_list.size() << endl;
						}
					}
					else {
						max_common_lcp = buffer[i].lcp;
					}
					top_lcp = buffer[i].lcp;

				}

				if(successive_color != current_color) {
					while(!stacked_list.empty()) {
						stackedLCPInterval wr_lcp_interval = stacked_list.back();
						stacked_list.pop_back();
						 cout << "D[" << wr_lcp_interval.pos <<
								"] = " << wr_lcp_interval.lcp + 1 << endl;
						 cout << stacked_list.size() << endl;
						fwrite(&wr_lcp_interval, sizeof(lcp_interval), 1, d_file);
						max_common_lcp = wr_lcp_interval.lcp;
					}
					current_color = successive_color;
				}

			} else {

				while(!stacked_list.empty()) {
					stacked_list.pop_front();
					// cout << stacked_list.size() << endl;
				}
				current_color = successive_color;
				max_common_lcp = 0;
				top_lcp = 0;

			}

		}
	}

	while ((num_read = fread(buffer, sizeof(buffer[0]), BUFFER_SIZE, f_ESA)) > 0)  {

		for(LetterNumber i = 0; i < num_read; ++i, ++k) {

			successive_color = mapColor(buffer[i].text);

			if(buffer[i].lcp != 0) {

				if(buffer[i].lcp > top_lcp) {

					lcp_interval = {
							k - 1,				// pos
							buffer[i].lcp,		// lcp value
					};
					stacked_list.push_front(lcp_interval);
					// cout << stacked_list.size() << endl;
					top_lcp = buffer[i].lcp;

				} else if (buffer[i].lcp < top_lcp) {

					LetterNumber ini_pos = i - 1;
					while(!stacked_list.empty() && buffer[i].lcp < top_lcp) {
						ini_pos = stacked_list.front().pos;
						stacked_list.pop_front();
						// cout << stacked_list.size() << endl;
						if(!stacked_list.empty())
							top_lcp = stacked_list.front().lcp;
						else
							top_lcp = 0;
					}

					if(buffer[i].lcp > max_common_lcp) {

						if(buffer[i].lcp > top_lcp) {

							lcp_interval = {
									ini_pos,			// pos
									buffer[i].lcp,	// lcp value
							};
							stacked_list.push_front(lcp_interval);
							// cout << stacked_list.size() << endl;
						}
					}
					else {
						max_common_lcp = buffer[i].lcp;
					}
					top_lcp = buffer[i].lcp;

				}

				if(successive_color != current_color) {
					while(!stacked_list.empty()) {
						stackedLCPInterval wr_lcp_interval = stacked_list.back();
						stacked_list.pop_back();
						// cout << "D[" << wr_lcp_interval.pos <<
						//		"] = " << wr_lcp_interval.lcp + 1 << endl;
						// cout << stacked_list.size() << endl;
						fwrite(&wr_lcp_interval, sizeof(lcp_interval), 1, d_file);
						max_common_lcp = wr_lcp_interval.lcp;
					}
					current_color = successive_color;
				}

			} else {

				while(!stacked_list.empty()) {
					stacked_list.pop_front();
					// cout << stacked_list.size() << endl;
				}
				current_color = successive_color;
				max_common_lcp = 0;
				top_lcp = 0;

			}

		}

	}

	fclose(f_ESA);
	fclose(d_file);

	cout << "D Computation End\n";
	cout << "Current memory: " << malloc_count_current() << " bytes" << endl;

	elapsed = (clock()-start)/(double)(CLOCKS_PER_SEC);
	elapsed_wc = difftime(time(NULL),start_wc);

	cout << "### D Computation Reporting ###\n";
	cout << "Elapsed time: " << elapsed << " secs\n";
	cout << "Wall Clock time: " << elapsed_wc << " secs\n";
	cout << "Peak memory: " << malloc_count_peak() << " bytes\n";
	cout << endl;

}

void StackedDGenerator::generateD() {

	std::string input_file_name = params->target_collection_file_name + C_GESAExt;
	FILE* f_ESA = fopen(input_file_name.c_str(), "rb");
	if(f_ESA == nullptr) {
		std::ostringstream err_message;
		err_message << "Couldn't open file " << input_file_name;
		Error::stopWithError(C_StackedDGenerator_ClassName, __func__, err_message.str());
	}
//	std::cout << "Read from EGSA File: " << input_file_name << std::endl;
	GESAReader gesa_reader(f_ESA);

	std::string d_file_name = params->output_file_name + C_DynBlockFileExt;
	FILE* d_file = fopen(d_file_name.c_str(), "wb");
	if(d_file == nullptr) {
		std::ostringstream err_message;
		err_message << "Couldn't open file " << d_file_name;
		Error::stopWithError(C_StackedDGenerator_ClassName, __func__, err_message.str());
	}
//	std::cout << "Write to D File: " << d_file_name << std::endl;
	LCPWriter d_writer(d_file);

	list<stackedLCPInterval> stacked_list;
	stackedLCPInterval lcp_interval, wr_lcp_interval;
	LetterNumber k = 0, last_pos = 0;
	SequenceLength top_lcp = 0, max_common_lcp = 0, max_stack_size = 0;
	bool current_color, successive_color;
	t_GSA gesa_struct;

	clock_t start = clock();
	time_t start_wc = time(NULL);
	double elapsed, elapsed_wc;

	cout << "D Computation Starting\n";
	cout << "Current memory: " << malloc_count_current() << " bytes" << endl;

	while(gesa_reader.readGESAStruct(gesa_struct) && gesa_struct.lcp == 0){

		// cout << "LCP: " << gesa_struct.lcp << "\tID: " << gesa_struct.text << endl;

		++k;
	}

	// cout << "LCP: " << gesa_struct.lcp << "\tID: " << gesa_struct.text << endl;

	lcp_interval = {
			k - 1,					// pos
			gesa_struct.lcp,		// lcp value
	};
	stacked_list.push_front(lcp_interval);
	if(max_stack_size < stacked_list.size())
		max_stack_size = stacked_list.size();
	// cout << stacked_list.size() << endl;
	top_lcp = gesa_struct.lcp;
	current_color = mapColor(gesa_struct.text);
	++k;

	while(gesa_reader.readGESAStruct(gesa_struct)) {

		// cout << "LCP: " << gesa_struct.lcp << "\tID: " << gesa_struct.text << endl;

		successive_color = mapColor(gesa_struct.text);

		if(gesa_struct.lcp != 0) {

			if(gesa_struct.lcp > top_lcp) {

				lcp_interval = {
						k - 1,					// pos
						gesa_struct.lcp,		// lcp value
				};
				stacked_list.push_front(lcp_interval);
				if(max_stack_size < stacked_list.size())
					max_stack_size = stacked_list.size();
				// cout << stacked_list.size() << endl;
				top_lcp = gesa_struct.lcp;

			} else if (gesa_struct.lcp < top_lcp) {

				LetterNumber ini_pos = k - 1;
				while(!stacked_list.empty() && gesa_struct.lcp < top_lcp) {
					ini_pos = stacked_list.front().pos;
					stacked_list.pop_front();
					// cout << stacked_list.size() << endl;
					if(!stacked_list.empty())
						top_lcp = stacked_list.front().lcp;
					else
						top_lcp = 0;
				}

				if(gesa_struct.lcp > max_common_lcp) {

					if(gesa_struct.lcp > top_lcp) {

						lcp_interval = {
								ini_pos,			// pos
								gesa_struct.lcp,	// lcp value
						};
						stacked_list.push_front(lcp_interval);
						if(max_stack_size < stacked_list.size())
							max_stack_size = stacked_list.size();
						// cout << stacked_list.size() << endl;
					}
				}
				else {
					max_common_lcp = gesa_struct.lcp;
				}
				top_lcp = gesa_struct.lcp;

			}

			if(successive_color != current_color) {
				while(!stacked_list.empty()) {
					wr_lcp_interval = stacked_list.back();
					stacked_list.pop_back();
//					cout << "D[" << wr_lcp_interval.pos <<
//							"] = " << wr_lcp_interval.lcp + 1 << endl;
//					 cout << stacked_list.size() << endl;
					d_writer.writeLCPPair(wr_lcp_interval, last_pos);
					last_pos = wr_lcp_interval.pos;
					max_common_lcp = wr_lcp_interval.lcp;
				}
				current_color = successive_color;
			}

		} else {

			while(!stacked_list.empty()) {
				stacked_list.pop_front();
				// cout << stacked_list.size() << endl;
			}
			current_color = successive_color;
			max_common_lcp = 0;
			top_lcp = 0;

		}
		++k;

	}
	d_writer.fillWithZeros(k - 1 - last_pos);
	d_writer.flushLCPValues();

	fclose(f_ESA);
	fclose(d_file);

	cout << "D Computation End\n";
	cout << "Current memory: " << malloc_count_current() << " bytes" << endl;

	elapsed = (clock()-start)/(double)(CLOCKS_PER_SEC);
	elapsed_wc = difftime(time(NULL),start_wc);

	cout << "### D Computation Reporting ###\n";
	cout << "Elapsed time: " << elapsed << " secs\n";
	cout << "Wall Clock time: " << elapsed_wc << " secs\n";
	cout << "Peak memory: " << malloc_count_peak() << " bytes\n";
	cout << "Max Block List Size: " << max_stack_size << endl;
	cout << endl;

}

} /* namespace multi_acs */
