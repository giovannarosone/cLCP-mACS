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
#include "GESAConverter.h"
#include "Types.h"
#include "Tools.h"

namespace multi_acs {

GESAConverter::GESAConverter() {
	// TODO Auto-generated constructor stub

}

GESAConverter::~GESAConverter() {
	// TODO Auto-generated destructor stub
}

void GESAConverter::extractFromGESA(const std::string &file_path) {

		std::string input_file_name = file_path + C_GESAExt;
		FILE* f_ESA;
		f_ESA = fopen(input_file_name.c_str(), "rb");
		if(f_ESA == nullptr) {
			std::ostringstream err_message;
			err_message << "Couldn't open file " << input_file_name;
			Error::stopWithError(C_GESAConverter_ClassName, __func__, err_message.str());
		}
		// std::cout << "Read EGSA File: " << input_file_name << std::endl;


		std::string output_file_name = file_path;
		FileName ebwt_file_name(output_file_name, C_BwtFileExt);
		FILE* ebwt_file = fopen(ebwt_file_name.c_str(), "wb");
		if(ebwt_file == nullptr) {
			std::ostringstream err_message;
			err_message << "Couldn't open file " << ebwt_file_name.str();
			Error::stopWithError(C_GESAConverter_ClassName, __func__, err_message.str());
		}
		// std::cout << "Write EBWT File: " << ebwt_file_name.str() << std::endl;

		FileName lcp_file_name(output_file_name, C_LcpFileExt);
		FILE* lcp_file = fopen(lcp_file_name.c_str(), "wb");
		if(lcp_file == nullptr) {
			std::ostringstream err_message;
			err_message << "Couldn't open file " << lcp_file_name.str();
			Error::stopWithError(C_GESAConverter_ClassName, __func__, err_message.str());
		}
		// std::cout << "Write LCP File: " << lcp_file_name.str() << std::endl;

		FileName id_file_name(output_file_name, C_IdFileExt);
		FILE* id_file = fopen(id_file_name.c_str(), "wb");
		if(id_file == nullptr) {
			std::ostringstream err_message;
			err_message << "Couldn't open file " << id_file_name.str();
			Error::stopWithError(C_GESAConverter_ClassName, __func__, err_message.str());
		}
		// std::cout << "Write Id File: " << id_file_name.str() << std::endl;

		t_GSA GSA;
		AlphabetSymbol c;
		SequenceNumber id;
		SequenceLength lcp;
		int num_read;
		while ((num_read = fread(&GSA, sizeof(t_GSA), 1, f_ESA)) > 0 )  {
			// cout << GSA.text << "\t" << GSA.suff << "\t" << GSA.lcp << "\t" << GSA.bwt << "\n";

			if (GSA.bwt == '\0')
				c = TERMINATE_CHAR;
			else
				c = GSA.bwt;
			fwrite(&c, sizeof(AlphabetSymbol), 1, ebwt_file);

			lcp = GSA.lcp;
			fwrite(&lcp, sizeof(SequenceLength), 1, lcp_file);

			id = GSA.text;
			fwrite(&id, sizeof(SequenceNumber), 1, id_file);

		}

		fclose(f_ESA);
		fclose(ebwt_file);
		fclose(lcp_file);
		fclose(id_file);

		std::cout << "Extraction of BWT/LCP/ID from " << input_file_name << " SUCCEEDED" << std::endl;
	}

} /* namespace multi_acs */


