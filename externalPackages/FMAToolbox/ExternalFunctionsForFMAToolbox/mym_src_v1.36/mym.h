/*  mYm is a Matlab interface to MySQL server that support BLOB object
*   Copyright (C) 2005 Swiss Federal Institute of technology (EPFL), Lausanne, CH
*
*   This program is free software; you can redistribute it and/or modify
*   it under the terms of the GNU General Public License as published by
*   the Free Software Foundation; either version 2 of the License, or
*   at your option) any later version.
*
*   This program is distributed in the hope that it will be useful,
*   but WITHOUT ANY WARRANTY; without even the implied warranty of
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*   GNU General Public License for more details.
*
*   You should have received a copy of the GNU General Public License
*   along with this program; if not, write to the Free Software
*   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301  USA
*
*   Author:                     Yannick Maret
*   e-mail:                     yannick.maret@a3.epfl.ch
*   old fashioned paper mail:   EPFL-STI-ITS-LTS1
*                               Yannick MARET
*                               ELD241
*                               Station 11
*                               CH-1015 Lausanne
*
*   Notice: some parts of this code (server connection, fancy print) is based  
*   on an original code by Robert Almgren (http://www.mmf.utoronto.ca/resrchres/mysql/).
*   The present code is under GPL license with the express agreement of Mr. Almgren.
*/

#ifndef MY_MAT_H
#define MY_MAT H

// some local defintion
#ifndef ulong
typedef unsigned long ulong;
#endif
#ifndef uchar
typedef unsigned char uchar;
#endif

#include <mex.h>  //  Definitions for Matlab API

#define ZLIB_WINAPI 
#include <zlib.h>

#include <math.h>

const bool debug = false;  //  turn on information messages

#if (defined(_WIN32)||defined(_WIN64))&&!defined(__WIN__)
#include <windows.h>
#include <winsock.h> // to overcome a bug in mysql.h
/* We use case-insensitive string comparison functions strcasecmp(), strncasecmp().
   These are a BSD addition and are also defined on Linux, but not on every OS, 
   in particular Windows. The two "inline" declarations below fix this problem. If
   you get errors on other platforms, move the declarations outside the WIN32 block */
inline int strcasecmp(const char *s1, const char *s2) { return strcmp(s1, s2); }
inline int strncasecmp(const char *s1, const char *s2, size_t n) { return strncmp(s1, s2, n); }
#endif
#include <mysql.h>  //  Definitions for MySQL client API

// DLL entry point
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
/*************** serializing functions ***************/
// matlab objects
char* serializeArray(ulong &rnBytes, const mxArray* rpArray, const char* rpArg, const bool);
char* serializeStruct(ulong &rnBytes, const mxArray* rpArray, const char* rpArg, const bool);
char* serializeCell(ulong &rnBytes, const mxArray* rpArray, const char* rpArg, const bool);
char* serializeScalar(ulong &rnBytes, const mxArray* rpArray, const char* rpArg, const bool);
char* serializeString(ulong &rnBytes, const mxArray* rpArray, const char* rpArg, const bool);
// generic object
char* serializeBinary(ulong &rnBytes, const mxArray* rpArray, const char* rpArg, const bool);
char* serializeFile(ulong &rnBytes, const mxArray* rpArray, const char* rpArg, const bool);
// deserilaizing functions
mxArray* deserialize(const char* rpSerial, const unsigned long rlength);
mxArray* deserializeArray(const char* rpSerial, const unsigned long rlength);
mxArray* deserializeStruct(const char* rpSerial, const unsigned long rlength);
mxArray* deserializeCell(const char* rpSerial, const unsigned long rlength);
// utility functioms
int file_length(FILE *f); // get the size of a file in byte
unsigned long min_mysql_escape(char* rpout, const char* rpin, const unsigned long nin);
/**********************************************************************
 *
 * hostport(s):  Given a host name s, possibly containing a port number
 *  separated by the port separation character (normally ':').
 * Modify the input string by putting a null at the end of
 * the host string, and return the integer of the port number.
 * Return zero in most special cases and error cases.
 * Modified string will not contain the port separation character.
 * Examples:  s = "myhost:2515" modifies s to "myhost" and returns 2515.
 *      s = "myhost"    leaves s unchanged and returns 0.
 *
 **********************************************************************/
const char portsep = ':';   //  separates host name from port number
static int hostport(char *s) {
	//   Look for first portsep in s; return 0 if null or can't find
	if (!s||!(s = strchr(s, portsep))) 
    return 0;
	//  If s points to portsep, then truncate and convert tail
	*s = 0;
	s+=1;
	return atoi(s);   // Returns zero in most special cases
}
/*********************************************************************/
//  Static variables that contain connection state
/*
 *  isopen gets set to true when we execute an "open"
 *  isopen gets set to false when either we execute a "close"
 *            or when a ping or status fails
 *   We do not set it to false when a normal query fails;
 *   this might be due to the server having died, but is much
 *   more likely to be caused by an incorrect query.
 */
class conninfo {
public:
	MYSQL *conn;   //  MySQL connection information structure
	bool isopen;   //  whether we believe that connection is open
	conninfo() { 
    conn = NULL; 
    isopen = false; 
  }
};
const int MAXCONN = 20;
static conninfo c[MAXCONN];   //  preserve state for MAXCONN different connections
// for GPL license
static bool runonce = false;
/*********************************************************************/
const char ID_MATLAB[] = "mYm";
const size_t LEN_ID_MATLAB = strlen(ID_MATLAB);
const char ID_ARRAY   = 'A';
const char ID_CELL    = 'C';
const char ID_STRUCT  = 'S';
// Placeholder related constants
const char PH_OPEN[]   = "{";  // placeholder openning symbols						
const char PH_CLOSE[]  = "}";  // placeholder closing symbols					
const char PH_BINARY   = 'B';
const char PH_FILE     = 'F';
const char PH_MATLAB   = 'M';
const char PH_STRING   = 'S';
// Preamble argument
const char PRE_NO_COMPRESSION = 'u';
// Compression
const ulong MIN_LEN_ZLIB = 1000;        // minimum number of byte for trying compressiom
const float MIN_CMP_RATE_ZLIB = 1.1f;   // minimum compression ratio
const char ZLIB_ID[] = "ZL123";
const size_t LEN_ZLIB_ID = strlen(ZLIB_ID);
typedef char* (*pfserial)(ulong &, const mxArray*, const char*, const bool);

static void getSerialFct(const char* rpt, const mxArray* rparg, pfserial& rpf, bool& rpec) {
	const unsigned n_dims = mxGetNumberOfDimensions(rparg);
	const int* p_dim = mxGetDimensions(rparg);
  bool no_compression = false;
  const char* pt = rpt;
  // first check for preamble
  if (*pt==PRE_NO_COMPRESSION) {
    no_compression = true;
    pt++;
  }
  if (*pt==PH_MATLAB) {
    // this placeholder results in a serialized version of array, cell or structure 
    rpec = true;
    if (mxIsNumeric(rparg)||mxIsChar(rparg))
      rpf = &serializeArray;
    else if (mxIsCell(rparg))
      rpf = &serializeCell;
    else if (mxIsStruct(rparg))
      rpf = &serializeStruct;
    else
      mexErrMsgTxt("Matlab placeholder only support array, structure, or cell");
  }
  else if (*pt==PH_BINARY||*pt==PH_FILE) {
    // this placeholder results in a binary dump of the corresponding data
    rpec = true;
    if (*pt==PH_BINARY)
      if (n_dims!=2||!(p_dim[0]==1||p_dim[1]==1)||!mxIsUint8(rparg))
				mexErrMsgTxt("Binary placeholders only accept UINT8 1-by-M or M-by-1 arrays!");
      else
        rpf = &serializeBinary;
    else
      if (n_dims!=2||!(p_dim[0]==1||p_dim[1]==1)||!mxIsChar(rparg))
        mexErrMsgTxt("String placeholders only accept CHAR 1-by-M or M-by-1 arrays!");
      else
        rpf = &serializeFile;
  }
  else if (*pt==PH_STRING) {
    // this placeholder results in a text version of the data
    rpec = false;
    rpf = &serializeString;
  }
  else
    mexErrMsgTxt("Unknow placeholders!");
  if (no_compression)
    rpec = false;
}
// entry point
mxArray* deserialize(const char* rpSerial, const unsigned long rlength) {
	mxArray* p_res = NULL;
	bool could_not_deserialize = true;
  bool used_compression = false;
  char* p_cmp = NULL;
  const char* p_serial = rpSerial;
  unsigned long length = rlength;
  if (p_serial==0) {
		// the row is empty: return an empty array
		p_res = mxCreateNumericArray(0, 0, mxCHAR_CLASS, mxREAL);
    return p_res;
  }
  if (strcmp(p_serial, ZLIB_ID)==0) {
    p_serial = p_serial+LEN_ZLIB_ID+1;
    // read the length in bytes
    int len;
    ulong lenLong;
    memcpy((char*)&len, p_serial, sizeof(int));
    p_serial = p_serial+sizeof(int);
    char* p_cmp = (char*)mxCalloc(len, sizeof(char));
    lenLong=(ulong)len;
    try {
      int res = uncompress((Bytef*)p_cmp, &lenLong, (const Bytef*)p_serial, length);
      if (res==Z_OK) {
        used_compression = true;
        p_serial = p_cmp;
        length = len;
      }
      else
        p_serial = rpSerial;
    }
    catch(...) {
      p_serial = rpSerial;
    }
  }
  if (strcmp(p_serial, ID_MATLAB)==0) {
    p_serial = p_serial+LEN_ID_MATLAB+1;
    try {
      could_not_deserialize = false;
      if (*p_serial==ID_ARRAY)
			  // the blob contains an array
			  p_res = deserializeArray(p_serial+1, length);
      else if (*p_serial==ID_STRUCT)
        // the blob contains a struct
        p_res = deserializeStruct(p_serial+1, length);
      else if (*p_serial==ID_CELL)
        // the blob contains a struct
        p_res = deserializeCell(p_serial+1, length);
      else
        // the blob contains an unknow object
        could_not_deserialize = true;
    }
    catch(...) {
      could_not_deserialize = true;
      p_serial = rpSerial;
    }
  }
	if (could_not_deserialize) {
		// we don't know what to do with the stuff in the database
		// we return a memory dump in UINT8 format
		int p_dim[2] = {length, 1};
		p_res = mxCreateNumericArray(2, p_dim, mxUINT8_CLASS, mxREAL);
		memcpy((char*)mxGetData(p_res), p_serial, length);
	}
  if (used_compression)
    mxFree(p_cmp);
	return p_res;
}
#endif // MY_MAT_H