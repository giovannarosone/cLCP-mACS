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
#ifndef COLLECTIONINFO_H_
#define COLLECTIONINFO_H_

#include "Types.h"
#include <string>

namespace multi_acs {

const char C_CollectionInfo_ClassName[]{"CollectionInfo"};
const std::string C_InfoExt{".info"};

class CollectionInfo {
public:
	CollectionInfo();
	CollectionInfo(const CollectionInfo &collection);
	CollectionInfo(const CollectionInfo &collection1, const CollectionInfo &collection2);
	CollectionInfo(const std::string &collection_file_name,
			const int input_format,
			const bool collect_info,
			const bool verbose);
	CollectionInfo(const std::string &collection_file_name);
	CollectionInfo& operator=(const CollectionInfo& collection);
	~CollectionInfo();

	LetterNumber size{0};
	std::map<AlphabetSymbol, LetterNumber> freq;
	std::map<SequenceNumber, SequenceLength> colors;

	AlphabetSize getAlphabetSize();
	SequenceNumber getSequenceNumber();
	SequenceLength getSequenceLength(const SequenceNumber id);
	void printCollectionInfo();
	void saveCollectionInfo();
	void loadCollectionInfo(std::string info_file_name);
	void loadCollectionLengths(std::string info_file_name);

	void join(const CollectionInfo &collection);

private:
	const std::string collection_file_name;
	const int file_format{0};
	const bool collect_info{true};
	const bool verbose{false};

	void checkForFiles();
	void checkForBCRFiles();
	void checkForGESAFiles();
	void collectSymbolsInfo(FILE* file);
	void collectColorsInfo(FILE* file);
	void collectSymbolsAndColorsInfo(FILE* file);
};

} /* namespace multi_acs */

#endif /* COLLECTIONINFO_H_ */
