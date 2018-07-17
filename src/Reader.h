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
#ifndef READER_H_
#define READER_H_

#include "Types.h"

namespace multi_acs {

class GESAReader {
public:
	GESAReader(FILE* gesa_file);
	~GESAReader();
	bool readGESAStruct(t_GSA &gesa_struct);

private:
	FILE* gesa_file{nullptr};
	t_GSA buffer[BUFFER_SIZE];
	LetterNumber struct_counter{BUFFER_SIZE};
	LetterNumber last_num_read{BUFFER_SIZE};
};

class EBWTReader {
public:
	EBWTReader(FILE* ebwt_file);
	~EBWTReader();
	bool readEBWTSymbol(AlphabetSymbol &symbol);

private:
	FILE* ebwt_file{nullptr};
	AlphabetSymbol buffer[BUFFER_SIZE]{'\0'};
	LetterNumber symbol_counter{BUFFER_SIZE};
	LetterNumber last_num_read{BUFFER_SIZE};
};

class IdReader {
public:
	IdReader(FILE* id_file);
	~IdReader();
	bool readSequenceId(SequenceNumber &id);

private:
	FILE* id_file{nullptr};
	SequenceNumber buffer[BUFFER_SIZE]{0};
	LetterNumber id_counter{BUFFER_SIZE};
	LetterNumber last_num_read{BUFFER_SIZE};
};

class LCPReader {
public:
	LCPReader(FILE* lcp_file);
	~LCPReader();
	bool readLCPValue(SequenceLength &value);

private:
	FILE* lcp_file{nullptr};
	SequenceLength buffer[BUFFER_SIZE]{0};
	LetterNumber lcp_counter{BUFFER_SIZE};
	LetterNumber last_num_read{BUFFER_SIZE};
};

class DReader {
public:
	DReader(FILE* d_file);
	~DReader();
	bool readDValue(SequenceLength &value);

private:
	FILE* d_file{nullptr};
	SequenceLength buffer[BUFFER_SIZE]{0};
	LetterNumber d_counter{BUFFER_SIZE};
	LetterNumber last_num_read{BUFFER_SIZE};
};

class ZReader {
public:
	ZReader(FILE* z_file);
	~ZReader();
	bool readZColor(SequenceNumber &color);

private:
	FILE* z_file{nullptr};
	SequenceNumber buffer[BUFFER_SIZE]{0};
	LetterNumber color_counter{BUFFER_SIZE};
	LetterNumber last_num_read{BUFFER_SIZE};
};

class BReader {
public:
	BReader(FILE* b_file);
	~BReader();
	bool readBlock(Block &block);

private:
	FILE* b_file{nullptr};
	Block buffer[BUFFER_SIZE];
	LetterNumber block_counter{BUFFER_SIZE};
	LetterNumber last_num_read{BUFFER_SIZE};
};

class IrrReader {
public:
	IrrReader(FILE* irr_file, const LetterNumber irr_file_size);
	~IrrReader();
	bool readIrrilevantBit(bool &irrilevant_bit);

private:
	FILE* irr_file{nullptr};
	const LetterNumber irr_file_size;
	LetterNumber global_bit_counter{0};
	BitGroup buffer[BUFFER_SIZE]{0};
	LetterNumber bit_counter{BUFFER_SIZE*(8*sizeof(BitGroup))};
	LetterNumber last_num_read{BUFFER_SIZE*(8*sizeof(BitGroup))};
};

} /* namespace multi_acs */

#endif /* READER_H_ */
