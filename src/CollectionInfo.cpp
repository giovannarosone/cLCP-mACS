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
#include "CollectionInfo.h"
#include "Tools.h"
#include <string>

using namespace std;

namespace multi_acs {

CollectionInfo::CollectionInfo() {

}

CollectionInfo::CollectionInfo(const std::string& collection_file_name,
		const int file_format,
		const bool verbose) :
	collection_file_name(collection_file_name),
	file_format(file_format),
	collect_info(false),
	verbose(verbose) {
	checkForFiles();
	loadCollectionInfo(collection_file_name);
}

CollectionInfo::CollectionInfo(const CollectionInfo& collection1,
		const CollectionInfo& collection2) {
	size = collection1.size;
	freq = collection1.freq;
	colors = collection1.colors;
	this->join(collection2);
}

CollectionInfo::CollectionInfo(const CollectionInfo& collection) :
	size(collection.size),
	freq(collection.freq),
	colors(collection.colors) { }

CollectionInfo& CollectionInfo::operator =(const CollectionInfo& collection) {
	size = collection.size;
	freq = collection.freq;
	colors = collection.colors;
	return *this;
}

CollectionInfo::CollectionInfo(const string &collection_file_name,
		const int file_format,
		const bool collect_info,
		const bool verbose) :
		collection_file_name(collection_file_name),
		file_format(file_format),
		collect_info(collect_info),
		verbose(verbose) {
	checkForFiles();
}

CollectionInfo::~CollectionInfo() {
}

AlphabetSize CollectionInfo::getAlphabetSize() {
	return freq.size();
}

SequenceNumber CollectionInfo::getSequenceNumber() {
	return colors.size();
}

SequenceLength CollectionInfo::getSequenceLength(const SequenceNumber id) {
	return colors[id];
}

void CollectionInfo::printCollectionInfo() {

	cout << "Collection Size (with separators): " << size << '\n';
	cout << "Collection Size (without separators): "
			<< size - colors.size() << '\n';
	cout << "Alphabet dimension: " << freq.size()
			<< " (" << freq.size() - 1 << " + " << TERMINATE_CHAR << ")\n";
	cout << "Symbols Frequency Distribution: ";
	for(auto it = freq.begin(); it != freq.end(); ++it) {
		if(it->first)
			cout << " [" << it->first << ':' << it->second << ']';
		else
			cout << " [" << TERMINATE_CHAR << ':' << it->second << ']';
	}
	cout << endl;
	cout << "Number of Sequences: " << colors.size() << '\n';
	cout << "Sequences Length Distribution\n";
	for(auto it = colors.begin(); it != colors.end(); ++it) {
		cout << "Seq " << it->first << ": " << it->second - 1 << '\n';
	}
	cout << endl;
}

void CollectionInfo::saveCollectionInfo() {

	FileName info_file_name(collection_file_name, C_InfoExt);
	FILE* info_file = fopen(info_file_name.c_str(), "w");
	if(info_file == nullptr) {
		std::ostringstream err_message;
		err_message << "Couldn't open file " << info_file;
		Error::stopWithError(C_CollectionInfo_ClassName, __func__, err_message.str());
	}

	fprintf(info_file, "%lu\n", size);
	fprintf(info_file, "#%lu\n", freq.size());
	for(auto it = freq.begin(); it != freq.end(); ++it) {
		fprintf(info_file, "%c\t%lu\n", it->first, it->second);
	}
	fprintf(info_file, "#%lu\n", colors.size());
	for(auto it = colors.begin(); it != colors.end(); ++it) {
		fprintf(info_file, "%lu\t%lu\n", it->first, it->second);
	}

	fclose(info_file);
	if(verbose) {
		cout << info_file_name.str() << " SAVED\n";
	}

}

void CollectionInfo::loadCollectionInfo(std::string file_name) {

	FileName info_file_name(file_name, C_InfoExt);
	FILE* info_file = fopen(info_file_name.c_str(), "r");
	if(info_file == nullptr) {
		std::ostringstream err_message;
		err_message << "Couldn't open file " << info_file;
		Error::stopWithError(C_CollectionInfo_ClassName, __func__, err_message.str());
	}

	fscanf(info_file, "%lu\n", &size);
	AlphabetSize alpha = 0;
	fscanf(info_file, "#%lu\n", &alpha);
	AlphabetSymbol symbol = '\0';
	LetterNumber frequency = 0;
	for(AlphabetSize i = 0; i < alpha; ++i) {
		fscanf(info_file, "%c\t%lu\n", &symbol, &frequency);
		freq.insert(pair<AlphabetSymbol, LetterNumber>(symbol, frequency));
	}
	SequenceNumber seqs_num = 0;
	fscanf(info_file, "#%lu\n", &seqs_num);
	SequenceNumber color = 0;
	SequenceLength length = 0;
	for(SequenceNumber i = 0; i < seqs_num; ++i) {
		fscanf(info_file, "%lu\t%lu\n", &color, &length);
		colors.insert(pair<SequenceNumber, SequenceLength>(color, length));
	}

	fclose(info_file);
	if(verbose) {
		cout << info_file_name.str() << " LOADED\n";
	}

}

void CollectionInfo::loadCollectionLengths(std::string info_file_name) {

		string len_file_name = info_file_name + ".lenSeqs.aux";
		FILE* len_file = fopen(len_file_name.c_str(), "rb");
		if (len_file == nullptr) {
			std::ostringstream err_message;
			err_message << "Couldn't open file " << len_file_name;
			Error::stopWithError(C_CollectionInfo_ClassName, __func__, err_message.str());
		}

		LetterNumber num_read = 0;
		SequenceNumber color = 0;
		vector<int_lcp> len_buffer(BUFFER_SIZE, 0);
		while((num_read = fread(&len_buffer[0], sizeof(int_lcp), len_buffer.size() , len_file)) > 0) {
			for(LetterNumber i = 0; i < num_read; ++i) {
				colors.insert(pair<SequenceNumber, SequenceLength>(color, len_buffer[i] + 1));
				++color;
				size += len_buffer[i] + 1;
			}
		}
		fclose(len_file);
}

void CollectionInfo::checkForFiles() {

	switch(file_format) {
		case 0: 	checkForBCRFiles();
					break;
		case 1:		checkForGESAFiles();
					break;
		default: 	ostringstream err_message;
					err_message << "Unrecognized file format";
					Error::stopWithError(C_CollectionInfo_ClassName, __func__, err_message.str());
	}

}

void CollectionInfo::checkForGESAFiles() {

	FileName gesa_file_name(collection_file_name, C_GESAExt);
	FILE* gesa_file = fopen(gesa_file_name.c_str(), "rb");
	if(gesa_file == nullptr) {
		std::ostringstream err_message;
		err_message << "Couldn't open file " << gesa_file;
		Error::stopWithError(C_CollectionInfo_ClassName, __func__, err_message.str());
	}
	if(collect_info) collectSymbolsAndColorsInfo(gesa_file);
	fclose(gesa_file);

	if(verbose) {
		cout << gesa_file_name.str() << " FOUND\n";
	}
}

void CollectionInfo::checkForBCRFiles() {

	FileName ebwt_file_name(collection_file_name, C_BwtFileExt);
	FILE* ebwt_file = fopen(ebwt_file_name.c_str(), "rb");
	if(ebwt_file == nullptr) {
		ostringstream err_message;
		err_message << "Couldn't open file " << ebwt_file_name.str();
		Error::stopWithError(C_CollectionInfo_ClassName, __func__, err_message.str());
	}
	if(collect_info) collectSymbolsInfo(ebwt_file);
	fclose(ebwt_file);

	if(verbose) {
		cout << ebwt_file_name.str() << " FOUND\n";
	}


	FileName lcp_file_name(collection_file_name, C_LcpFileExt);
	FILE* lcp_file = fopen(lcp_file_name.c_str(), "rb");
	if(lcp_file == nullptr) {
		ostringstream err_message;
		err_message << "Couldn't open file " << lcp_file_name.str();
		Error::stopWithError(C_CollectionInfo_ClassName, __func__, err_message.str());
	}
	fclose(lcp_file);

	if(verbose) {
		cout << lcp_file_name.str() << " FOUND\n";
	}


	FileName id_file_name(collection_file_name, C_IdFileExt);
	FILE* id_file = fopen(id_file_name.c_str(), "rb");
	if(id_file == nullptr) {
		ostringstream err_message;
		err_message << "Couldn't open file " << id_file_name.str();
		Error::stopWithError(C_CollectionInfo_ClassName, __func__, err_message.str());
	}
	if(collect_info) collectColorsInfo(id_file);
	fclose(id_file);

	if(verbose) {
		cout << id_file_name.str() << " FOUND\n";
	}

}

void CollectionInfo::collectSymbolsInfo(FILE* file) {

	AlphabetSymbol buffer[BUFFER_SIZE];
	LetterNumber num_read = 0;
	while((num_read = fread(buffer, sizeof(AlphabetSymbol), BUFFER_SIZE, file)) > 0) {
		size += num_read;
		for(LetterNumber i = 0; i < num_read; ++i) {
			if(freq.count(buffer[i]) > 0) {
				++freq[buffer[i]];
			} else {
				freq.insert(pair<AlphabetSymbol, LetterNumber>(buffer[i], 1));
			}
		}
	}

}

void CollectionInfo::collectColorsInfo(FILE* file) {

	SequenceNumber buffer[BUFFER_SIZE];
	LetterNumber num_read = 0;
	while((num_read = fread(buffer, sizeof(SequenceNumber), BUFFER_SIZE, file)) > 0) {
		for(LetterNumber i = 0; i < num_read; ++i) {
			if(colors.count(buffer[i]) > 0) {
				++colors[buffer[i]];
			} else {
				colors.insert(pair<SequenceNumber, SequenceLength>(buffer[i], 1));
			}
		}
	}

}

void CollectionInfo::collectSymbolsAndColorsInfo(FILE* file) {

	t_GSA buffer[BUFFER_SIZE];
	LetterNumber num_read = 0;
	while((num_read = fread(buffer, sizeof(t_GSA), BUFFER_SIZE, file)) > 0) {
		size += num_read;
		for(LetterNumber i = 0; i < num_read; ++i) {
			if(freq.count(buffer[i].bwt) > 0) {
				++freq[buffer[i].bwt];
			} else {
				freq.insert(pair<AlphabetSymbol, LetterNumber>(buffer[i].bwt, 1));
			}
			if(colors.count(buffer[i].text) > 0) {
				++colors[buffer[i].text];
			} else {
				colors.insert(pair<SequenceNumber, SequenceLength>(buffer[i].text, 1));
			}
		}
	}

}

void CollectionInfo::join(const CollectionInfo& collection) {

	size += collection.size;
	for(const pair<AlphabetSymbol, LetterNumber> &freq_pair : collection.freq) {
		if(freq.count(freq_pair.first) > 0) {
			freq[freq_pair.first] += freq_pair.second;
		} else {
			freq.insert(pair<AlphabetSymbol, LetterNumber>(freq_pair.first, freq_pair.second));
		}
	}
	SequenceLength sequence_number = colors.size();
	for(const pair<SequenceNumber, SequenceLength> &colors_pair : collection.colors) {
		if(colors.count(colors_pair.first) > 0) {
			colors.insert(pair<SequenceNumber, SequenceLength>(colors_pair.first + sequence_number, colors_pair.second));
		} else {
			freq.insert(pair<SequenceNumber, SequenceLength>(colors_pair.first, colors_pair.second));
		}
	}

}

} /* namespace multi_acs */
