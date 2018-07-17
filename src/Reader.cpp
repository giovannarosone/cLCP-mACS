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
#include "Reader.h"

namespace multi_acs {

GESAReader::GESAReader(FILE* gesa_file) :
	gesa_file(gesa_file) { }

GESAReader::~GESAReader() {

}

bool GESAReader::readGESAStruct(t_GSA &gesa_struct) {
	LetterNumber num_read = 0;
	if(struct_counter == last_num_read) {
		num_read = fread(buffer, sizeof(t_GSA), BUFFER_SIZE, gesa_file);
		if(num_read == 0) {
			return false;
		} else {
			last_num_read = num_read;
		}
		struct_counter = 0;
	}
	gesa_struct = buffer[struct_counter++];
	return true;
}


EBWTReader::EBWTReader(FILE* ebwt_file) :
	ebwt_file(ebwt_file) { }

EBWTReader::~EBWTReader() {

}

bool EBWTReader::readEBWTSymbol(AlphabetSymbol& symbol) {
	LetterNumber num_read = 0;
	if(symbol_counter == last_num_read) {
		num_read = fread(buffer, sizeof(AlphabetSymbol), BUFFER_SIZE, ebwt_file);
		if(num_read == 0) {
			return false;
		} else {
			last_num_read = num_read;
		}
		symbol_counter = 0;
	}
	symbol = buffer[symbol_counter++];
	return true;
}

IdReader::IdReader(FILE* id_file) :
	id_file(id_file) { }

IdReader::~IdReader() {

}

bool IdReader::readSequenceId(SequenceNumber& id) {
	LetterNumber num_read = 0;
	if(id_counter == last_num_read) {
		num_read = fread(buffer, sizeof(SequenceNumber), BUFFER_SIZE, id_file);
		if(num_read == 0) {
			return false;
		} else {
			last_num_read = num_read;
		}
		id_counter = 0;
	}
	id = buffer[id_counter++];
	return true;
}

LCPReader::LCPReader(FILE* lcp_file) :
	lcp_file(lcp_file) { }

LCPReader::~LCPReader() {

}

bool LCPReader::readLCPValue(SequenceLength& value) {
	LetterNumber num_read = 0;
	if(lcp_counter == last_num_read) {
		num_read = fread(buffer, sizeof(SequenceLength), BUFFER_SIZE, lcp_file);
		if(num_read == 0) {
			return false;
		} else {
			last_num_read = num_read;
		}
		lcp_counter = 0;
	}
	value = buffer[lcp_counter++];
	return true;
}

DReader::DReader(FILE* d_file) :
	d_file(d_file) { }

DReader::~DReader() {

}

bool DReader::readDValue(SequenceLength& value) {
	LetterNumber num_read = 0;
	if(d_counter == last_num_read) {
		num_read = fread(buffer, sizeof(SequenceLength), BUFFER_SIZE, d_file);
		if(num_read == 0) {
			return false;
		} else {
			last_num_read = num_read;
		}
		d_counter = 0;
	}
	value = buffer[d_counter++];
	return true;
}

ZReader::ZReader(FILE* z_file) :
	z_file(z_file) { }

ZReader::~ZReader() {

}

bool ZReader::readZColor(SequenceNumber& color) {
	LetterNumber num_read = 0;
	if(color_counter == last_num_read) {
		num_read = fread(buffer, sizeof(SequenceNumber), BUFFER_SIZE, z_file);
		if(num_read == 0) {
			return false;
		} else {
			last_num_read = num_read;
		}
		color_counter = 0;
	}
	color = buffer[color_counter++];
	return true;
}

BReader::BReader(FILE* b_file) :
	b_file(b_file) { }

BReader::~BReader() {

}

bool BReader::readBlock(Block& block) {
	LetterNumber num_read = 0;
	if(block_counter == last_num_read) {
		num_read = fread(buffer, sizeof(Block), BUFFER_SIZE, b_file);
		if(num_read == 0) {
			return false;
		} else {
			last_num_read = num_read;
		}
		block_counter = 0;
	}
	block = buffer[block_counter++];
	return true;
}

IrrReader::IrrReader(FILE* irr_file, const LetterNumber irr_file_size) :
	irr_file(irr_file),
	irr_file_size(irr_file_size){ }

IrrReader::~IrrReader() {

}

bool IrrReader::readIrrilevantBit(bool &irrilevant_bit) {
	LetterNumber num_read = 0;
	if(bit_counter == last_num_read) {
		num_read = fread(buffer, sizeof(BitGroup), BUFFER_SIZE, irr_file);
		if(num_read == 0) {
			return false;
		} else {
			last_num_read = num_read*8;
		}
		bit_counter = 0;
	}
	if(global_bit_counter < irr_file_size) {
		int byte_pos = bit_counter/(8*sizeof(BitGroup));
		int shift_in_byte = bit_counter%(8*sizeof(BitGroup));
		irrilevant_bit = buffer[byte_pos] & (1 << shift_in_byte);
		++bit_counter;
		++global_bit_counter;
		return true;
	}
	return false;
}

} /* namespace multi_acs */

