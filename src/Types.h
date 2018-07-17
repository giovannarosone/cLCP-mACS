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
#ifndef TYPES_H_
#define TYPES_H_

#include <iostream>
#include <vector>
#include <map>

namespace multi_acs {

// Type to represent: Group of bit to read/write in a bitfile
typedef unsigned char BitGroup;

// Type to represent: Symbol of the sequences alphabet (Typically {A, C, G T})
typedef unsigned char AlphabetSymbol;

// Type to represent: Dimension of the sequences alphabet
typedef AlphabetSymbol AlphabetSize;

/* TODO Handle correct size of the type used.
 * - Is it better to provide a macro variable to specify the types dimension at compile time
 */

// Type to represent: Length of a sequence / LCP value
typedef uint32_t SequenceLength;
// typedef uint16_t SequenceLength;

// Type to represent: Number of sequences
/* USE: uint32_t - below 4.294.967.296 sequences
 *		uint64_t - otherwise
 */
typedef uint32_t SequenceNumber;

// Type to represent: Number of characters in EBWT / SA value (position)
/* USE: uint32_t - below 4.294.967.296 characters
 *		uint64_t - otherwise
 */
typedef uint64_t LetterNumber;

// Type to represent: Value of distance between sequences
typedef double DistanceValue;

// Type to represent: Maximum memory usage for computation
typedef uint64_t AllocableMemory;

// Element of the generalized suffix array (GSA)
struct GSAElement {
	// Position in the sequence
	SequenceLength pos;
	// Number of the sequence
	SequenceNumber seq;
};

struct Block {
	SequenceLength round;
	LetterNumber size;
	bool is_monochrome;
	SequenceNumber color;
};

struct stackedLCPInterval {
	LetterNumber pos;
	SequenceLength lcp;
};

// eGSA Types

typedef unsigned int int_text;
typedef unsigned int int_suff;
typedef unsigned int int_lcp;
typedef unsigned char int8;

#pragma pack(1)
typedef struct{

	int_text	text;
	int_suff	suff;
	int_lcp 	lcp;
	int8		bwt;

} t_GSA;

const std::string C_GESAExt{".gesa"};

#define BUFFER_SIZE 10000

#define PLACEHOLDER_CHAR '#'
#define TERMINATE_CHAR '$'
// #define TERMINATE_CHAR '\0'

#define FIRST_COLOR 0
#define SECOND_COLOR 1

const std::string C_BwtFileExt{".bwt"};
const std::string C_LcpFileExt{".lcp"};
const std::string C_SaFileExt{".sa"};
const std::string C_IdFileExt{".id"};
const std::string C_BlockFileExt{".b"};
const std::string C_IrrBlockBitfileExt{".ir"};
const std::string C_DynBlockFileExt{".d"};
const std::string C_ColorFileExt{".z"};
const std::string C_CLcpFileExt{".clcp"};
const std::string C_PartialCLcpFileExt{".xclcp"};
const std::string C_DistanceFileExt{".acs"};

const AlphabetSymbol C_MaxAlphabetSize{static_cast<AlphabetSize>(-1)};
const SequenceNumber C_MaxSequenceNumber{static_cast<SequenceNumber>(-1)};
const SequenceLength C_MaxSequenceLength{static_cast<SequenceLength>(-1)};
const LetterNumber C_MaxLetterNumber{static_cast<LetterNumber>(-1)};



struct UnresolvedBlock {
	bool is_open{false};
	bool is_monochrome{true};
	SequenceNumber color{C_MaxSequenceNumber};
	LetterNumber start{C_MaxLetterNumber};
	LetterNumber end{C_MaxLetterNumber};

	void reset() {
		is_open = false;
		is_monochrome = true;
		color = C_MaxSequenceNumber;
		start = C_MaxLetterNumber;
		end = C_MaxLetterNumber;
	}
};

} /* namespace multi_acs */

#endif /* TYPES_H_ */
