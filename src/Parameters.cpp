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
#include "Parameters.h"
#include "Types.h"
#include <iostream>
#include <string>

using namespace std;

namespace multi_acs {

// Gap ------------------------------------------

GapParameters::GapParameters(const bool verbose,
		const string first_file_name,
		const string second_file_name,
		const int file_format,
		const string output_file_name) :
				verbose(verbose),
				first_collection_file_name(first_file_name),
				second_collection_file_name(second_file_name),
				file_format(file_format),
				output_file_name(output_file_name) { }

GapParameters::~GapParameters() {
	// TODO Auto-generated destructor stub
}

void GapParameters::printParameters() {

	cout << "+-------------------- Resume Options --- ---------------+\n";
	cout << "First Collection file path: "
			<< first_collection_file_name << '\n';
	cout << "First Collection file path: "
			<< second_collection_file_name << '\n';
	cout << "Output file path: "
			<< output_file_name << '\n';
	if(verbose)
		cout << "***WITH VERBOSE REPORTING***\n";
	cout << "+-------------------------------------------------------+" << endl;
}

// MultiACS -------------------------------------

MultiACSParameters::MultiACSParameters(const bool verbose,
		const string reference_sequence_file_name,
		const string target_collection_file_name,
		const int file_format,
		const SequenceNumber reference_color,
		const string output_file_name,
		const AllocableMemory memory_amount) :
				verbose(verbose),
				reference_sequence_file_name(reference_sequence_file_name),
				target_collection_file_name(target_collection_file_name),
				file_format(file_format),
				reference_color(reference_color),
				output_file_name(output_file_name),
				memory_amount(memory_amount) { }

MultiACSParameters::~MultiACSParameters() {
	// TODO Auto-generated destructor stub
}

void MultiACSParameters::printParameters() {

	cout << "----- OPTIONS RESUME -----\n";
	cout << "Reference Sequence file path: "
			<< reference_sequence_file_name << '\n';
	cout << "Target Collection file path: "
			<< target_collection_file_name << '\n';
	cout << "Reference Sequence Color: "
			<< reference_color << '\n';
	cout << "Output file path: "
			<< output_file_name << '\n';
	cout << "Max Usable Memory in ACS Computing: " << memory_amount << " Byte\n";
	if(verbose)
		cout << "***WITH VERBOSE REPORTING***\n";
	cout << "--- END OPTIONS RESUME ---" << endl;
}

} /* namespace multi_acs */
