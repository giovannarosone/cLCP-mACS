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
#ifndef WRITER_H_
#define WRITER_H_

#include "Types.h"

namespace multi_acs {

const char C_Writer_ClassName[]{"Writer"};

class InplaceZSegmentWriter {
public:
	InplaceZSegmentWriter(int z_segment_file);
	~InplaceZSegmentWriter();
	void writeColor(SequenceNumber color);
	void skipColors(LetterNumber color_to_skip);
	void flushColors();

private:
	int z_segment_file{-1};
	SequenceNumber buffer[BUFFER_SIZE]{'\0'};
	LetterNumber color_counter{0};
	LetterNumber z_file_offset{0};
	LetterNumber last_read{0};
};

class ZSegmentWriter {
public:
	ZSegmentWriter(FILE* z_segment_file);
	~ZSegmentWriter();
	void writeColor(SequenceNumber color);
	void flushColors();

private:
	FILE* z_segment_file{nullptr};
	SequenceNumber buffer[BUFFER_SIZE]{0};
	LetterNumber color_counter{0};
};

class IdWriter {
public:
	IdWriter(FILE* id_file);
	~IdWriter();
	void writeSequenceId(SequenceNumber id);
	void flushIds();

private:
	FILE* id_file{nullptr};
	SequenceNumber buffer[BUFFER_SIZE]{0};
	LetterNumber id_counter{0};
};

class LCPWriter {
public:
	LCPWriter(FILE* lcp_file);
	~LCPWriter();
	void writeLCPValue(SequenceLength lcp_value);
	void writeLCPPair(stackedLCPInterval lcp_interval, LetterNumber last_pos);
	void fillWithZeros(LetterNumber zeros_num);
	void flushLCPValues();

private:
	FILE* lcp_file{nullptr};
	SequenceLength buffer[BUFFER_SIZE];
	LetterNumber lcp_counter{0};
};

class InplaceBSegmentWriter {
public:
	InplaceBSegmentWriter(int b_file);
	~InplaceBSegmentWriter();
	SequenceLength writeBValue(SequenceLength b_value);
	void skipBValues(LetterNumber b_to_skip);
	void flushBValues();

private:
	int b_file{-1};
	SequenceLength buffer[BUFFER_SIZE]{0};
	LetterNumber b_counter{0};
	LetterNumber b_file_offset{0};
	LetterNumber last_read{0};
};

class IrrSegmentWriter {
public:
	IrrSegmentWriter(FILE* irr_segment_file);
	~IrrSegmentWriter();
	void writeIrrilevantBit(bool irrilevant_bit);
	void flushBits();

private:
	FILE* irr_segment_file{nullptr};
	LetterNumber irr_segment_size{0};
	BitGroup buffer[BUFFER_SIZE]{0};
	LetterNumber bit_counter{0};
};

class InplaceIrrSegmentWriter {
public:
	InplaceIrrSegmentWriter(int irr_segment_file);
	~InplaceIrrSegmentWriter();
	void writeIrrilevantBit(bool irrilevant_bit);
	void skipBits(LetterNumber bit_to_skip);
	void flushBits();

private:
	int irr_segment_file{-1};
	LetterNumber irr_segment_size{0};
	BitGroup buffer[BUFFER_SIZE]{0};
	LetterNumber bit_counter{0};
	LetterNumber irr_file_offset{0};
	LetterNumber last_read{0};
};

class EBWTWriter {
public:
	EBWTWriter(FILE* ebwt_file);
	~EBWTWriter();
	void writeSymbol(AlphabetSymbol symbol);
	void flushSymbols();

private:
	FILE* ebwt_file{nullptr};
	AlphabetSymbol buffer[BUFFER_SIZE]{'\0'};
	LetterNumber symbol_counter{0};
};

class CLCPWriter {
public:
	CLCPWriter(FILE* clcp_file, SequenceNumber seq_num, SequenceLength buffer_size = 0);
	~CLCPWriter();
	void writeCLCPValues(SequenceLength *clcp_values);
	void flushCLCPValues();

private:
	FILE* clcp_file{nullptr};
	SequenceNumber seq_num{1};
	SequenceLength buffer_size{BUFFER_SIZE};
	SequenceLength **buffer{nullptr};
	LetterNumber clcp_counter{0};
};


} /* namespace multi_acs */

#endif /* WRITER_H_ */
