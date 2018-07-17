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
#include "Writer.h"
#include "Tools.h"
#include <unistd.h>
#include <sstream>

using namespace std;

namespace multi_acs {

// Z ----------------------------------------------------------

InplaceZSegmentWriter::InplaceZSegmentWriter(int z_segment_file) :
		z_segment_file(z_segment_file) {
}

InplaceZSegmentWriter::~InplaceZSegmentWriter() {

}

void InplaceZSegmentWriter::writeColor(SequenceNumber color) {
	if(color_counter == last_read) {
		LetterNumber num_write = pwrite(z_segment_file, buffer, sizeof(SequenceNumber)*last_read, z_file_offset);
		z_file_offset += num_write;
		color_counter = 0;
		last_read = pread(z_segment_file, buffer, sizeof(SequenceNumber)*BUFFER_SIZE, z_file_offset);
		last_read /= sizeof(SequenceNumber);
	}
	buffer[color_counter] = color;
	++color_counter;
}

void InplaceZSegmentWriter::skipColors(LetterNumber color_to_skip) {
	if(color_counter + color_to_skip > last_read) {
		LetterNumber remaining = last_read - color_counter;
		color_counter += remaining;
		color_to_skip -= remaining;
		flushColors();
		z_file_offset += sizeof(SequenceNumber)*color_to_skip;
		last_read = pread(z_segment_file, buffer, sizeof(SequenceNumber)*BUFFER_SIZE, z_file_offset);
		last_read /= sizeof(SequenceNumber);
	}
	else {
		color_counter += color_to_skip;
	}
}

void InplaceZSegmentWriter::flushColors() {
	if(color_counter > 0) {
		LetterNumber num_write = pwrite(z_segment_file, buffer, sizeof(SequenceNumber)*color_counter, z_file_offset);
		z_file_offset += num_write;
		color_counter = 0;
	}
}



ZSegmentWriter::ZSegmentWriter(FILE* z_segment_file) :
		z_segment_file(z_segment_file) { }

ZSegmentWriter::~ZSegmentWriter() {

}

void ZSegmentWriter::writeColor(SequenceNumber color) {
	if(color_counter == BUFFER_SIZE) {
		fwrite(buffer, sizeof(SequenceNumber), BUFFER_SIZE, z_segment_file);
		color_counter = 0;
	}
	buffer[color_counter] = color;
	++color_counter;
}

void ZSegmentWriter::flushColors() {
	if(color_counter) {
		fwrite(buffer, sizeof(SequenceNumber), color_counter, z_segment_file);
		color_counter = 0;
	}
}

// Id ---------------------------------------------------------

IdWriter::IdWriter(FILE* id_file) :
		id_file(id_file) { }

IdWriter::~IdWriter() {

}

void IdWriter::writeSequenceId(SequenceNumber id) {
	if(id_counter == BUFFER_SIZE) {
		fwrite(buffer, sizeof(SequenceNumber), BUFFER_SIZE, id_file);
		id_counter = 0;
	}
	buffer[id_counter] = id;
	++id_counter;
}

void IdWriter::flushIds() {
	if(id_counter) {
		fwrite(buffer, sizeof(SequenceNumber), id_counter, id_file);
		id_counter = 0;
	}
}


// Irr --------------------------------------------------------

IrrSegmentWriter::IrrSegmentWriter(FILE* irr_segment_file) :
		irr_segment_file(irr_segment_file) { }

IrrSegmentWriter::~IrrSegmentWriter() {

}

void IrrSegmentWriter::writeIrrilevantBit(bool irrilevant_bit) {
	if(bit_counter == BUFFER_SIZE*8*sizeof(BitGroup)) {
		fwrite(buffer, sizeof(BitGroup), BUFFER_SIZE, irr_segment_file);
		bit_counter = 0;
	}
	int byte_pos = bit_counter/(8*sizeof(BitGroup));
	int shift_in_byte = bit_counter%(8*sizeof(BitGroup));
	buffer[byte_pos] |= ((irrilevant_bit?1:0) << shift_in_byte);
	++bit_counter;
}

void IrrSegmentWriter::flushBits() {
	if(bit_counter) {
		int byte_num = bit_counter/(8*sizeof(BitGroup)) + 1;
		fwrite(buffer, sizeof(BitGroup), byte_num, irr_segment_file);
		bit_counter = 0;
	}
}

InplaceIrrSegmentWriter::InplaceIrrSegmentWriter(int irr_segment_file) :
		irr_segment_file(irr_segment_file) {
}

InplaceIrrSegmentWriter::~InplaceIrrSegmentWriter() {

}

void InplaceIrrSegmentWriter::writeIrrilevantBit(bool irrilevant_bit) {
	if(bit_counter == last_read*sizeof(BitGroup)*8) {
		LetterNumber num_write = pwrite(irr_segment_file, buffer, sizeof(BitGroup)*last_read, irr_file_offset);
		irr_file_offset += num_write;
		bit_counter = 0;
		last_read = pread(irr_segment_file, buffer, sizeof(BitGroup)*BUFFER_SIZE, irr_file_offset);
		last_read /= sizeof(BitGroup);
	}
	int byte_pos = bit_counter/(8*sizeof(BitGroup));
	int shift_in_byte = bit_counter%(8*sizeof(BitGroup));
	buffer[byte_pos] |= ((irrilevant_bit?1:0) << shift_in_byte);
	++bit_counter;
	// if(irrilevant_bit) cout << 1 ;
}

void InplaceIrrSegmentWriter::skipBits(LetterNumber bit_to_skip) {
	if(bit_counter + bit_to_skip > last_read*8*sizeof(BitGroup)) {
		int remaining_bit = (sizeof(BitGroup)*8 - bit_counter%8)%8;
		bit_counter += remaining_bit;
		bit_to_skip -= remaining_bit;
		flushBits();
		int byte_num = bit_to_skip/(8*sizeof(BitGroup));
		int excedent_bits = bit_to_skip%(8*sizeof(BitGroup));
		irr_file_offset += sizeof(BitGroup)*byte_num;
		last_read = pread(irr_segment_file, buffer, sizeof(BitGroup)*BUFFER_SIZE, irr_file_offset);
		last_read /= sizeof(BitGroup);
		bit_counter += excedent_bits;
	}
	else {
		bit_counter += bit_to_skip;
	}
}

void InplaceIrrSegmentWriter::flushBits() {
	if(bit_counter > 0) {
		int byte_num = ceil((double)bit_counter/(8*sizeof(BitGroup)));
		LetterNumber num_write = pwrite(irr_segment_file, buffer, sizeof(BitGroup)*byte_num, irr_file_offset);
		irr_file_offset += num_write;
		bit_counter = 0;
	}
}

// EBWT -------------------------------------------------------

EBWTWriter::EBWTWriter(FILE* ebwt_file) :
		ebwt_file(ebwt_file) { }

EBWTWriter::~EBWTWriter() {

}

void EBWTWriter::writeSymbol(AlphabetSymbol symbol) {
	if(symbol_counter == BUFFER_SIZE) {
		fwrite(buffer, sizeof(AlphabetSymbol), BUFFER_SIZE, ebwt_file);
		symbol_counter = 0;
	}
	buffer[symbol_counter] = symbol;
	++symbol_counter;
}

void EBWTWriter::flushSymbols() {
	fwrite(buffer, sizeof(AlphabetSymbol), symbol_counter, ebwt_file);
	symbol_counter = 0;
}

// LCP/B ------------------------------------------------------

LCPWriter::LCPWriter(FILE* lcp_file) :
		lcp_file(lcp_file) { }

LCPWriter::~LCPWriter() {

}

void LCPWriter::writeLCPValue(SequenceLength lcp_value) {
	if(lcp_counter == BUFFER_SIZE) {
		fwrite(buffer, sizeof(SequenceLength), BUFFER_SIZE, lcp_file);
		lcp_counter = 0;
	}
	buffer[lcp_counter] = lcp_value;
	++lcp_counter;
}

void LCPWriter::flushLCPValues() {
	if(lcp_counter) {
		fwrite(buffer, sizeof(SequenceLength), lcp_counter, lcp_file);
		lcp_counter = 0;
	}
}

void LCPWriter::fillWithZeros(LetterNumber zeros_num) {
	if(lcp_counter + zeros_num > BUFFER_SIZE) {
		LetterNumber remaining = BUFFER_SIZE - lcp_counter;
		zeros_num -= remaining;
		for(LetterNumber i = 0; i < remaining; ++i)
			buffer[lcp_counter++] = 0;
		fwrite(buffer, sizeof(SequenceLength), lcp_counter, lcp_file);
		lcp_counter = 0;
		while(zeros_num > 0) {
			if(zeros_num < BUFFER_SIZE) {
				for(LetterNumber i = 0; i < zeros_num; ++i)
					buffer[lcp_counter++] = 0;
				zeros_num = 0;
			}
			else {
				for(LetterNumber i = 0; i < BUFFER_SIZE; ++i)
					buffer[i] = 0;
				fwrite(buffer, sizeof(SequenceLength), BUFFER_SIZE, lcp_file);
				zeros_num -= BUFFER_SIZE;
			}
		}
	}
	else {
		for(LetterNumber i = 0; i < zeros_num; ++i)
			buffer[lcp_counter++] = 0;
	}
}

void LCPWriter::writeLCPPair(stackedLCPInterval lcp_interval, LetterNumber last_pos) {

	LetterNumber to_write = lcp_interval.pos - last_pos + (last_pos ? 0 : 1);
	if(lcp_counter + to_write > BUFFER_SIZE) {
		LetterNumber remaining = BUFFER_SIZE - lcp_counter;
		to_write -= remaining;
		for(LetterNumber i = 0; i < remaining; ++i)
			buffer[lcp_counter++] = 0;
		fwrite(buffer, sizeof(SequenceLength), lcp_counter, lcp_file);
		lcp_counter = 0;
		while(to_write > 0) {
			if(to_write < BUFFER_SIZE) {
				for(LetterNumber i = 0; i < to_write - 1; ++i)
					buffer[lcp_counter++] = 0;
				buffer[lcp_counter++] = lcp_interval.lcp + 1;
				to_write = 0;
			}
			else {
				for(LetterNumber i = 0; i < BUFFER_SIZE; ++i)
					buffer[i] = 0;
				fwrite(buffer, sizeof(SequenceLength), BUFFER_SIZE, lcp_file);
				to_write -= BUFFER_SIZE;
			}
		}

	}
	else {
		for(LetterNumber i = 0; i < to_write - 1; ++i)
			buffer[lcp_counter++] = 0;
		buffer[lcp_counter++] = lcp_interval.lcp + 1;
	}

}

InplaceBSegmentWriter::InplaceBSegmentWriter(int b_file) :
		b_file(b_file) {
}

InplaceBSegmentWriter::~InplaceBSegmentWriter() {

}

SequenceLength InplaceBSegmentWriter::writeBValue(SequenceLength b_value) {
	if(b_counter == last_read) {
		LetterNumber num_write = pwrite(b_file, buffer, sizeof(SequenceLength)*last_read, b_file_offset);
		b_file_offset += num_write;
		b_counter = 0;
		last_read = pread(b_file, buffer, sizeof(SequenceLength)*BUFFER_SIZE, b_file_offset);
		last_read /= sizeof(SequenceLength);
	}
	SequenceLength old_b_value = buffer[b_counter];
	if(old_b_value == 0)
		buffer[b_counter] = b_value;
	++b_counter;
	return old_b_value;
}

void InplaceBSegmentWriter::skipBValues(LetterNumber b_to_skip) {
	if(b_counter + b_to_skip > last_read) {
		LetterNumber remaining = last_read - b_counter;
		b_counter += remaining;
		b_to_skip -= remaining;
		flushBValues();
		b_file_offset += b_to_skip*sizeof(SequenceLength);
		last_read = pread(b_file, buffer, sizeof(SequenceLength)*BUFFER_SIZE, b_file_offset);
		last_read /= sizeof(SequenceLength);
	}
	else {
		b_counter += b_to_skip;
	}
}

void InplaceBSegmentWriter::flushBValues() {
	if(b_counter > 0) {
		LetterNumber num_write = pwrite(b_file, buffer, sizeof(SequenceLength)*b_counter, b_file_offset);
		b_file_offset += num_write;
		b_counter = 0;
	}
}

CLCPWriter::CLCPWriter(FILE* clcp_file, SequenceNumber seq_num, SequenceLength buffer_size) :
		clcp_file(clcp_file),
		seq_num(seq_num),
		buffer_size(buffer_size) {

	buffer = new SequenceLength*[buffer_size]{nullptr};
	for(SequenceLength i = 0; i < buffer_size; ++i) {
		buffer[i] = new SequenceLength[seq_num]{0};
	}

}

CLCPWriter::~CLCPWriter() {
	for(SequenceLength i = 0; i < buffer_size; ++i) {
		delete[] buffer[i];
	}
	delete[] buffer;
}

void CLCPWriter::writeCLCPValues(SequenceLength* lcp_values) {
	if(clcp_counter == buffer_size) {
		fwrite(buffer, sizeof(SequenceLength), buffer_size*seq_num, clcp_file);
		clcp_counter = 0;
	}
	for(SequenceNumber i = 0; i < seq_num; ++i)
		buffer[clcp_counter][i] = lcp_values[i];
	++clcp_counter;
}

void CLCPWriter::flushCLCPValues() {
	if(clcp_counter) {
		fwrite(buffer, sizeof(SequenceLength), clcp_counter*seq_num, clcp_file);
		clcp_counter = 0;
	}
}


}


