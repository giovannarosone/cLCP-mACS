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
#ifndef TOOLS_H_
#define TOOLS_H_

#include "Types.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <algorithm>

namespace multi_acs {

class FileName {
public:
	FileName( const std::string &str )
		: m_str( str ) { }

	~FileName() { }

	// EBWT/LCP/GSA/Z/B/D segments file name
	FileName(const int i, const std::string &ext) {
		std::ostringstream fn;
		fn << i << ext;
		m_str = fn.str();
	}

	// EBWT/LCP/GSA/Z/B/D file name
	FileName(const char *file_name, const std::string &ext) {
		std::ostringstream fn;
		fn << file_name << ext;
		m_str = fn.str();
	}

	FileName(const std::string &file_name, const std::string &ext) {
		m_str = file_name + ext;
	}

	// New EBWT/LCP/GSA/Z/B/D segments file name
	FileName(const std::string &file_name, const int i, const std::string &ext) {
		std::ostringstream fn;
		fn << file_name << i << ext;
		m_str = fn.str();
	}

	std::string str() const {
		return m_str;
	}

	const char* c_str() const {
		return m_str.c_str();
	}

	static std::string removeExtension(const std::string &file_name) {
		return file_name.substr(0, file_name.find_last_of("."));
	}

	static std::string removePath(const std::string &file_name) {
		return file_name.substr(file_name.find_last_of("/") + 1, file_name.length());
	}

	static std::string getPath(const std::string &file_name) {
		size_t i = file_name.find_last_of("/");
		if(i != std::string::npos)
			return file_name.substr(0, file_name.find_last_of("/"));
		else
			return ".";
	}

private:
	std::string m_str;
};

class ErrorMessage {
public:
	ErrorMessage(const std::string error_message) :
	m_error_message(error_message) { }
private:
	const std::string m_error_message;
};

class Error {
public:
	static void stopWithError(const char* class_name, const char* funcion_name,
			const std::string &message) {
		std::ostringstream stream;
		stream << class_name << "::" << funcion_name << " ERROR: " << message;
		std::cerr << stream.str() << std::endl;
		exit(EXIT_FAILURE);
	}

	static void continueWithWarning(const char* class_name, const char* funcion_name,
			const std::string &message) {
		std::ostringstream stream;
		stream << class_name << "::" << funcion_name << " WARNING: " << message;
		std::cerr << stream.str() << std::endl;
	}
};

//bool compare_files(const std::string& filename1, const std::string& filename2)
//{
//    std::ifstream file1(filename1);
//    std::ifstream file2(filename2);
//
//    std::istreambuf_iterator<char> begin1(file1);
//    std::istreambuf_iterator<char> begin2(file2);
//
//    std::istreambuf_iterator<char> end;
//
//    return std::equal(begin1, end, begin2);
//}


} /* namespace multi_acs */

#endif /* TOOLS_H_ */
