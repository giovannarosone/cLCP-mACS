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
#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include "Types.h"
#include <string>

namespace multi_acs {


class GapParameters {
public:
	GapParameters(const bool verbose,
			const std::string first_file_name,
			const std::string second_file_name,
			const int file_format,
			const std::string output_file_name);
	~GapParameters();
	const bool verbose;
	const std::string first_collection_file_name;
	const std::string second_collection_file_name;
	const int file_format;
	const std::string output_file_name;
	void printParameters();
};

class MultiACSParameters {
public:
	MultiACSParameters(const bool verbose,
			const std::string reference_sequence_file_name,
			const std::string target_collection_file_name,
			const int file_format,
			const SequenceNumber reference_color,
			const std::string output_file_name,
			const AllocableMemory memory_amount);
	~MultiACSParameters();
	const bool verbose;
	const std::string reference_sequence_file_name;
	const std::string target_collection_file_name;
	const int file_format;
	const SequenceNumber reference_color;
	const std::string output_file_name;
	const AllocableMemory memory_amount;
	void printParameters();
};

} /* namespace multi_acs */

#endif /* PARAMETERS_H_ */
