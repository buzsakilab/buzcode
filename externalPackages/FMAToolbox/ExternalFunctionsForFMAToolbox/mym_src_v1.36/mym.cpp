/*  mYm is a Matlab interface to MySQL server that support BLOB object
 *   
 *   
 *
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <algorithm> // for max
#include "mym.h"

/**********************************************************************
 *typestr(s):  Readable translation of MySQL field type specifier
 *       as listed in   mysql_com.h
 **********************************************************************/
static const char*typestr(enum_field_types t) {
    switch(t) {
        //  These are considered numeric by IS_NUM() macro
        case FIELD_TYPE_DECIMAL:
            return "decimal";
        case FIELD_TYPE_TINY:
            return "tiny";
        case FIELD_TYPE_SHORT:
            return "short";
        case FIELD_TYPE_LONG:
            return "long";
        case FIELD_TYPE_FLOAT:
            return "float";
        case FIELD_TYPE_DOUBLE:
            return "double";
        case FIELD_TYPE_NULL:
            return "null";
        case FIELD_TYPE_LONGLONG:
            return "longlong";
        case FIELD_TYPE_INT24:
            return "int24";
        case FIELD_TYPE_YEAR:
            return "year";
        case FIELD_TYPE_TIMESTAMP:
            return "timestamp";
            //  These are not considered numeric by IS_NUM()
        case FIELD_TYPE_DATE:
            return "date";
        case FIELD_TYPE_TIME:
            return "time";
        case FIELD_TYPE_DATETIME:
            return "datetime";
        case FIELD_TYPE_NEWDATE:
            return "newdate";   // not in manual
        case FIELD_TYPE_ENUM:
            return "enum";
        case FIELD_TYPE_SET:
            return "set";
        case FIELD_TYPE_TINY_BLOB:
            return "tiny_blob";
        case FIELD_TYPE_MEDIUM_BLOB:
            return "medium_blob";
        case FIELD_TYPE_LONG_BLOB:
            return "long_blob";
        case FIELD_TYPE_BLOB:
            return "blob";
        case FIELD_TYPE_VAR_STRING:
            return "var_string";
        case FIELD_TYPE_STRING:
            return "string";
        case FIELD_TYPE_GEOMETRY:
            return "geometry";   // not in manual 4.0
            default:
                return "unknown";
    }
}
/**********************************************************************
 *fancyprint():  Print a nice display of a query result
 *  We assume the whole output set is already stored in memory,
 *  as from mysql_store_result(), just so that we can get the
 *  number of rows in case we need to clip the printing.
 *   In any case, we make only one pass through the data.
 *  If the number of rows in the result is greater than NROWMAX,
 *  then we print only the first NHEAD and the last NTAIL.
 *  NROWMAX must be greater than NHEAD+NTAIL, normally at least
 *  2 greater to allow the the extra information
 *  lines printed when we clip (ellipses and total lines).
 *   Display null elements as empty
 **********************************************************************/
const char*contstr = "...";   //  representation of missing rows
const int contstrlen = 3;     //  length of above string
const int NROWMAX = 50;       //  max number of rows to print w/o clipping
const int NHEAD = 10;         //  number of rows to print at top if we clip
const int NTAIL = 10;         //  number of rows to print at end if we clip
static void fancyprint(MYSQL_RES*res) {
    const ulong nrow = (ulong)mysql_num_rows(res);
    const ulong nfield = mysql_num_fields(res);
    bool clip = (nrow > NROWMAX);
    MYSQL_FIELD*f = mysql_fetch_fields(res);
    /************************************************************************/
    //  Determine column widths, and format strings for printing
    //  Find the longest entry in each column header,
    //  and over the rows, using MySQL's max_length
    ulong*len = (ulong*) mxMalloc(nfield*sizeof(ulong));
    for (ulong j = 0; j<nfield; j++) {
        len[j] = (int)strlen(f[j].name);
        if (f[j].max_length > len[j])
            len[j] = f[j].max_length;
    }
    //  Compare to the continuation string length if we will clip
    if (clip) {
        for (ulong j = 0; j<nfield; j++)
            if (contstrlen > len[j])
                len[j] = contstrlen;
    }
    //  Construct the format specification for printing the strings
    char**fmt = (char**) mxMalloc(nfield*sizeof(char*));
    for (ulong j = 0; j<nfield; j++) {
        fmt[j] = (char*) mxCalloc(10, sizeof(char));
        sprintf(fmt[j], "  %%-%ds ", len[j]);
    }
    /************************************************************************/
    //  Now print the actual data
    mexPrintf("\n");
    //  Column headers
    for (ulong j = 0; j<nfield; j++)
        mexPrintf(fmt[j], f[j].name);
    mexPrintf("\n");
    //  Fancy underlines
    for (ulong j = 0; j<nfield; j++) {
        mexPrintf(" +");
        for (ulong k = 0; k<len[j]; k++)
            mexPrintf("-");
        mexPrintf("+");
    }
    mexPrintf("\n");
    //  Print the table entries
    if (nrow<=0)
        mexPrintf("(zero rows in result set)\n");
    else {
        if (!clip) {
            //  print the entire table
            mysql_data_seek(res, 0);
            for (ulong i = 0; i<nrow; i++) {
                MYSQL_ROW row = mysql_fetch_row(res);
                if (!row) {
                    mexPrintf("Printing full table data from row %d\n", i+1);
                    mexErrMsgTxt("Internal error:  Failed to get a row");
                }
                for (ulong j = 0; j<nfield; j++)
                    mexPrintf(fmt[j], (row[j] ? row[j] : ""));
                    mexPrintf("\n"); }
        }
        else {
            //  print half at beginning, half at end
            mysql_data_seek(res, 0);
            for (int i = 0; i<NHEAD; i++) {
                MYSQL_ROW row = mysql_fetch_row(res);
                if (!row) {
                    mexPrintf("Printing head table data from row %d\n", i+1);
                    mexErrMsgTxt("Internal error:  Failed to get a row");
                }
                for (ulong j = 0; j<nfield; j++)
                    mexPrintf(fmt[j], (row[j] ? row[j] : ""));
                    mexPrintf("\n");
            }
            for (ulong j = 0; j<nfield; j++)
                mexPrintf(fmt[j], contstr);
            mexPrintf("\n");
            mysql_data_seek(res, nrow-NTAIL);
            for (int i = 0; i<NTAIL; i++) {
                MYSQL_ROW row = mysql_fetch_row(res);
                if (!row) {
                    mexPrintf("Printing tail table data from row %d", nrow-NTAIL+i+1);
                    mexErrMsgTxt("Internal error:  Failed to get a row");
                }
                for (ulong j = 0; j<nfield; j++)
                    mexPrintf(fmt[j], (row[j] ? row[j] : ""));
                    mexPrintf("\n");
            }
            mexPrintf("(%d rows total)\n", nrow);
        }
    }
    mexPrintf("\n");
    // These should be automatically freed when we return to Matlab,
    //  but just in case ...
    mxFree(len);
    mxFree(fmt);
}
/**********************************************************************
 *field2num():  Convert field in string format to double number
 **********************************************************************/
const double NaN = mxGetNaN();				//  Matlab NaN for null values
static bool can_convert(enum_field_types t) {
    return (IS_NUM(t)&&
    (t!=FIELD_TYPE_DATE)&&
    (t!=FIELD_TYPE_TIME)&&
    (t!=FIELD_TYPE_DATETIME)&&
    (t!=FIELD_TYPE_TIMESTAMP));
}
static double field2num(const char*s, enum_field_types t) {
    if (!s) return NaN;  // MySQL null -- nothing there
    if (IS_NUM(t)) {
        double val = NaN;
        char* charend = 0;
        val = std::strtod(s, &charend);
        if (charend==s) {
            //if (sscanf(s, "%lf", &val)!=1) {
            mexPrintf("Unreadable value \"%s\" of type %s\n", s, typestr(t));
            return NaN;
        }
        return val;
    }
    else {
        mexPrintf("Tried to convert \"%s\" of type %s to numeric\n", s, typestr(t));
        mexErrMsgTxt("Internal inconsistency");
    }
    return 0;
}
/**********************************************************************
 *getstring():   Extract string from a Matlab array
 *   (Space allocated by mxCalloc() should be freed by Matlab
 *  when control returns out of the MEX-function.)
 *  This is based on an original by Kimmo Uutela
 **********************************************************************/
static char* getstring(const mxArray*a) {
    int llen = mxGetM(a)*mxGetN(a)*sizeof(mxChar)+1;
    char*c = (char*) mxCalloc(llen, sizeof(char));
    if (mxGetString(a, c, llen))
        mexErrMsgTxt("Can\'t copy string in getstring()");
    return c;
}
/**********************************************************************
 *mysql():  Execute the actual action
 * Which action we perform is based on the first input argument,
 * which must be present and must be a character string:
 *   'open', 'close', 'use', 'status', or a legitimate MySQL query.
 *  This version does not permit binary characters in query string,
 * since the query is converted to a C null-terminated string.
 *  If no output argument is given, then information is displayed.
 * If an output argument is given, then we operate silently, and
 * return status information.
 **********************************************************************/
void mexFunction(int nlhs, mxArray*plhs[], int nrhs, const mxArray*prhs[]) {
    int cid = 0;   // ID number of the current connection
    int jarg = 0;  // Number of first string arg: becomes 1 if id is specified
    /*********************************************************************/
    // show GPL license
    if (!runonce) {
        runonce = true;
        mexPrintf("Matlab<->MySQL connector v1.36\n");
    }
    // Parse the first argument to see if it is a specific id number
    if (nrhs!=0&&mxIsNumeric(prhs[0]))	{
        if (mxGetM(prhs[0])!=1||mxGetN(prhs[0])!=1) {
            mexPrintf("Usage:  %s([id], command, [ host, user, password ])\n", mexFunctionName());
            mexPrintf("First argument is array %d x %d\n", mxGetM(prhs[0]), mxGetN(prhs[0]));
            mexErrMsgTxt("Invalid connection ID");
        }
        double xid = *mxGetPr(prhs[0]);
        cid = int(xid);
        if (double(cid)!=xid||cid<-1||cid>=MAXCONN) {
            mexPrintf("Usage:  %s([id], command, [ host, user, password ])\n", mexFunctionName());
            mexPrintf("id = %g -- Must be integer between 0 and %d\n", xid, MAXCONN-1);
            mexErrMsgTxt("Invalid connection ID");
        }
        jarg++;
    }
    if (debug)
        mexPrintf("cid = %d  jarg = %d\n", cid, jarg);
    /*********************************************************************/
    //  Parse the result based on the first argument
    enum querytype { OPEN, CLOSE, CLOSE_ALL, USE, STATUS, CMD } q;
    char*query = NULL;
    if (nrhs<=jarg)
        q = STATUS;
    else{
        if (!mxIsChar(prhs[jarg]))
            mexErrMsgTxt("The command string is missing!");
        query = getstring(prhs[jarg]);
        if (!strcasecmp(query, "open"))
            q = OPEN;
        else if (!strcasecmp(query, "close"))
            q = CLOSE;
        else if (!strcasecmp(query, "closeall"))
            q = CLOSE_ALL;
        else if (!strcasecmp(query, "use")||!strncasecmp(query, "use ", 4))
            q = USE;
        else if (!strcasecmp(query, "status"))
            q = STATUS;
        else q = CMD;
    }
    //  Check that the arguments are all character strings
    if (q!=CMD)
        for (int j = jarg; j<nrhs; j++) {
            if (!mxIsChar(prhs[j])) {
                mexPrintf("Usage:  %s([id], command, [ host, user, password ])\n", mexFunctionName());
                mexErrMsgTxt("All args must be strings, except id");
            }
        }
    if (debug){
        if (q==OPEN)
            mexPrintf("q = OPEN\n");
        if (q==CLOSE)
            mexPrintf("q = CLOSE\n");
        if (q==USE)
            mexPrintf("q = USE\n");
        if (q==STATUS)
            mexPrintf("q = STATUS\n");
        if (q==CMD)
            mexPrintf("q = CMD\n");
    }
    // found the correct cid, if automatic scan is requested
    if (cid==-1) {
        if (q!=OPEN)
            mexErrMsgTxt("automatic free connection scan, only when issuing a mYm open command!");
        cid = 1;
        while (cid<MAXCONN&&c[cid].isopen)
            cid++;
        if (cid==MAXCONN)
            mexErrMsgTxt("no free connection ID");
    }
    // Shorthand notation now that we know which id we have
    typedef MYSQL*mp;
    mp &conn = c[cid].conn;
    bool &isopen = c[cid].isopen;
    if (q==OPEN) {
        // OPEN A NEW CONNECTION
        ////////////////////////
        //  Close connection if it is open
        if (isopen) {
            mysql_close(conn);
            isopen = false;
            conn = NULL;
        }
        //  Extract information from input arguments
        char*host = NULL;
        if (nrhs>=jarg+2)
            host = getstring(prhs[jarg+1]);
        char*user = NULL;
        if (nrhs>=jarg+3)
            user = getstring(prhs[jarg+2]);
        char*pass = NULL;
        if (nrhs>=jarg+4)
            pass = getstring(prhs[jarg+3]);
        int port = hostport(host);  // returns zero if there is no port
        if (nlhs<1) {
            mexPrintf("Connecting to  host = %s", (host) ? host : "localhost");
            if (port)
                mexPrintf("  port = %d", port);
            if (user)
                mexPrintf("  user = %s", user);
            if (pass)
                mexPrintf("  password = %s", pass);
            mexPrintf("\n");
        }
        //  Establish and test the connection
        //  If this fails, then conn is still set, but isopen stays false
        if (!(conn = mysql_init(conn)))
            mexErrMsgTxt("Couldn\'t initialize MySQL connection object");
        //my_bool  my_true = true;
        //mysql_options(conn, MYSQL_OPT_RECONNECT, &my_true);
        if (!mysql_real_connect(conn, host, user, pass, NULL, port, NULL, CLIENT_MULTI_STATEMENTS))
            mexErrMsgTxt(mysql_error(conn));
        const char*c = mysql_stat(conn);
        if (c) {
            if (nlhs<1)
                mexPrintf("%s\n", c);
        }
        else
            mexErrMsgTxt(mysql_error(conn));
        isopen = true;
        
        //  Now we are OK -- if he wants output, give him the cid
        if (nlhs>=1) {
            if (!(plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL)))
                mexErrMsgTxt("Unable to create matrix for output");
            *mxGetPr(plhs[0]) = cid;
        }
    }
    else if (q==CLOSE) {
        // CLOSE CONNECTION
        ////////////////////////
        if (isopen) {
            // only if open
            mysql_close(conn);
            isopen = false;
            conn = NULL;
        }
    }
    else if (q==CLOSE_ALL) {
        // CLOSE ALL CONNECTIONS
        ////////////////////////
        for (int i = 0; i<MAXCONN; i++)
            if (c[i].isopen) {
                mysql_close(c[i].conn);
                c[i].conn = NULL;
                c[i].isopen = false;
            }
    }
    else if (q==USE) {
        // SPECIFY DATABASE
        ////////////////////////
        if (!isopen)
            mexErrMsgTxt("Not connected");
        if (mysql_ping(conn)) {
            isopen = false;
            mexErrMsgTxt(mysql_error(conn));
        }
        char*db = NULL;
        if (!strcasecmp(query, "use")) {
            if (nrhs>=2)
                db = getstring(prhs[jarg+1]);
            else
                mexErrMsgTxt("Must specify a database to use");
        }
        else if (!strncasecmp(query, "use ", 4)) {
            db = query+4;
            while (*db==' ' ||*db=='\t')
                db++;
        }
        else
            mexErrMsgTxt("How did we get here?  Internal logic error!");
        if (mysql_select_db(conn, db))
            mexErrMsgTxt(mysql_error(conn));
        if (nlhs<1)
            mexPrintf("Current database is \"%s\"\n", db);
        else {
            // he wants some output, give him a 1
            if (!(plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL)))
                mexErrMsgTxt("Unable to create matrix for output");
            *mxGetPr(plhs[0]) = 1.;
        }
    }
    else if (q==STATUS)	{
        // GET CONNECTION STATUS
        ////////////////////////
        if (nlhs<1) {
            //  He just wants a report
            if (jarg>0) {
                //  mysql(cid, 'status')  Give status of the specified id
                if (!isopen) {
                    mexPrintf("%2d: Not connected\n", cid);
                    return;
                }
                if (mysql_ping(conn)) {
                    isopen = false;
                    mexErrMsgTxt(mysql_error(conn));
                }
                mexPrintf("%2d:  %-30s   Server version %s\n", cid, mysql_get_host_info(conn), mysql_get_server_info(conn));
            }
            else {
                //  mysql('status') with no specified connection
                int nconn = 0;
                for (int j = 0; j<MAXCONN; j++)
                    if (c[j].isopen)
                        nconn++;
                if (debug)
                    mexPrintf("%d connections open\n", nconn);
                if (nconn==0) {
                    mexPrintf("No connections open\n");
                    return;
                }
                if (nconn==1&&c[0].isopen) {
                    // Only connection number zero is open
                    // Give simple report with no connection id #
                    if (mysql_ping(conn)) {
                        isopen = false;
                        mexErrMsgTxt(mysql_error(conn));
                    }
                    mexPrintf("Connected to %s Server version %s Client %s\n", mysql_get_host_info(conn), mysql_get_server_info(conn), mysql_get_client_info());
                }
                else {
                    // More than one connection is open
                    // Give a detailed report of all open connections
                    for (int j = 0; j<MAXCONN; j++) {
                        if (c[j].isopen) {
                            if (mysql_ping(c[j].conn)) {
                                c[j].isopen = false;
                                mexPrintf("%2d:  %s\n", mysql_error(c[j].conn));
                                continue;
                            }
                            mexPrintf("%2d: %-30s Server version %s\n", j, mysql_get_host_info(c[j].conn), mysql_get_server_info(c[j].conn));
                        }
                    }
                }
            }
        }
        else {
            //  He wants a return value for this connection
            if (!(plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL)))
                mexErrMsgTxt("Unable to create matrix for output");
            double *pr = mxGetPr(plhs[0]);
            *pr = 0.;
            if (!isopen)
                *pr = 1.;
            //return;
            else if (mysql_ping(conn)) {
                isopen = false;
                *pr = 2.;
                //return;
            }
            else if (!mysql_stat(conn)) {
                isopen = false;
                *pr = 3.;
                //return;
            }
        }
    }
    else if (q==CMD) {
        // ISSUE A MYSQL COMMAND
        ////////////////////////
        //  Check that we have a valid connection
        if (!isopen)
            mexErrMsgTxt("Not connected");
        if (mysql_ping(conn)) {
            isopen = false;
            mexErrMsgTxt(mysql_error(conn));
        }
        //******************PLACEHOLDER PROCESSING******************
        // global placeholders variables and constant
        const unsigned nex = nrhs-jarg-1;	// expected number of placeholders
        unsigned nac = 0;									// actual number
        ulong lq = (ulong)strlen(query);	// original query length, needed at some point
        if (nex) {
            // local placeholders variables and constant
            char** po = 0;								// pointer to placeholders openning symbols
            char** pc = 0;								// pointer to placeholders closing symbols
            bool*  pec = 0;               // pointer to 'enable compression field'
            char** pa = 0;								// pointer to placeholder in-string argument
            unsigned* ps = 0;							// pointer to placeholders size (in bytes)
            pfserial* pf = 0;							// pointer to serialization function
            // LOOK FOR THE PLACEHOLDERS
            po = (char**)mxCalloc(nex+1, sizeof(char*));
            pc = (char**)mxCalloc(nex+1, sizeof(char*));
            if (po[nac++] = strstr(query, PH_OPEN))
                while (po[nac-1]&&nac<=nex) {
                    pc[nac-1] = strstr(po[nac-1]+1, PH_CLOSE);
                    if (pc[nac-1]==0)
                        mexErrMsgTxt("Placeholders are not correctly closed!");
                    po[nac] = strstr(pc[nac-1]+1, PH_OPEN);
                    nac++;
                }
            nac--;
            // some error checks
            if (nac!=nex)
                mexErrMsgTxt("The number of placeholders differs from that of additional arguments!");
            // now we have the correct number of placeholders
            // read the types and in-string arguments
            ps = (unsigned*)mxCalloc(nac, sizeof(unsigned));
            pa = (char**)mxCalloc(nac, sizeof(char*));
            pec = (bool*)mxCalloc(nac, sizeof(bool));
            pf = (pfserial*)mxCalloc(nac, sizeof(pfserial));
            for (unsigned i = 0; i<nac; i++) {
                // placeholder function+control ph type+control arg type
                getSerialFct(po[i]+strlen(PH_OPEN), prhs[jarg+i+1], pf[i], pec[i]);
                // placeholder size in bytes
                ps[i] = (unsigned)(pc[i]-po[i]+strlen(PH_OPEN));
                // placeholder in-string argument
                pa[i] = po[i]+strlen(PH_OPEN)+1;
                *(pc[i]) = 0;
                // we can add a termination character in the query,
                // it will be processed anyway and we already have the original query length
            }
            // BUILD THE NEW QUERY
            // serialize
            char** pd = (char**)mxCalloc(nac, sizeof(char*));
            ulong* plen = (ulong*)mxCalloc(nac, sizeof(ulong));
            uLongf cmp_buf_len = 10000;
            Bytef* pcmp = (Bytef*)mxCalloc(cmp_buf_len, sizeof(Bytef));
            ulong nb = 0;
            for (unsigned i = 0; i<nac; i++) {
                // serialize individual field
                pd[i] = pf[i](plen[i], prhs[jarg+i+1], pa[i], true);
                if (pec[i]&&plen[i]>MIN_LEN_ZLIB) {
                    // compress field
                    const uLongf max_len = compressBound(plen[i]);
                    if (max_len>cmp_buf_len){
                        pcmp = (Bytef*)mxRealloc(pcmp, max_len);
                        cmp_buf_len = max_len;
                    }
                    uLongf len = cmp_buf_len;
                    if (compress(pcmp, &len, (Bytef*)pd[i], plen[i])!=Z_OK)
                        mexErrMsgTxt("Compression failed");
                    const float cmp_rate = plen[i]/(LEN_ZLIB_ID+1.f+sizeof(int)+len);
                    if (cmp_rate>MIN_CMP_RATE_ZLIB) {
                        const int len_old = (int)plen[i];
                        plen[i] = (ulong)LEN_ZLIB_ID+1+sizeof(int)+len;
                        pd[i] = (char*)mxRealloc(pd[i], plen[i]);
                        memcpy(pd[i], (const char*)ZLIB_ID, LEN_ZLIB_ID);
                        *(pd[i]+LEN_ZLIB_ID) = 0;
                        memcpy(pd[i]+LEN_ZLIB_ID+1, (const char*)&len_old, sizeof(int));
                        memcpy(pd[i]+LEN_ZLIB_ID+1+sizeof(int), pcmp, len);
                    }
                }
                nb += plen[i];
            }
            // create new query
            char* nq = (char*)mxCalloc(2*nb+lq+1, sizeof(char));  // new query
            char* pnq = nq; // running pointer to new query
            const char* poq = query; // running pointer to old query
            for (unsigned i = 0; i<nac; i++) {
                memcpy(pnq, poq, po[i]-poq);
                pnq += po[i]-poq;
                poq = po[i]+ps[i];
                pnq += mysql_real_escape_string(conn, pnq, pd[i], plen[i]);
                mxFree(pd[i]);
            }
            memcpy(pnq, poq, lq-(poq-query)+1);
            // replace the old query by the new one
            pnq += lq-(poq-query)+1;
            mxFree(query);
            query = nq;
            lq = (ulong)(pnq-nq);
            mxFree(po);
            mxFree(pc);
            mxFree(pec);
            mxFree(ps);
            mxFree(pa);
            mxFree(pd);
            mxFree(plen);
            if (pcmp!=NULL)
                mxFree(pcmp);
        }
        else {
            // check that no placeholders are present in the query
            char*p_tmp_o = strstr(query, PH_OPEN);
            char*p_tmp_c = strstr(query, PH_CLOSE);
            if (p_tmp_o||p_tmp_c)
                mexErrMsgTxt("The query contains placeholders, but no additional arguments!");
        }
        //  Execute the query (data stays on server)
        if (nac!=0) {
            if (mysql_real_query(conn, query, lq))
                mexErrMsgTxt(mysql_error(conn));
        }
        else if (mysql_query(conn, query))
            mexErrMsgTxt(mysql_error(conn));
        //  Download the data from server into our memory
        //  We need to be careful to deallocate res before returning.
        //  Matlab's allocation routines return instantly if there is not
        //  enough free space, without giving us time to dealloc res.
        //  This is a potential memory leak but I don't see how to fix it.
        MYSQL_RES*res = mysql_store_result(conn);
        //  As recommended in Paul DuBois' MySQL book (New Riders, 1999):
        //  A NULL result set after the query can indicate either
        //  (1) the query was an INSERT, DELETE, REPLACE, or UPDATE, that
        //    affect rows in the table but do not return a result set; or
        //  (2) an error, if the query was a SELECT, SHOW, or EXPLAIN
        //    that should return a result set but didn't.
        //  Distinguish between the two by checking mysql_field_count()
        //  We return in either case, either correctly or with an error
        if (!res) {
            if (!mysql_field_count(conn)) {
                ulong nrows = (ulong)mysql_affected_rows(conn);
                if (nlhs<1)
                    return;
                else {
                    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
                    if (plhs[0]==0)
                        mexErrMsgTxt("Unable to create numeric matrix for output");
                    *mxGetPr(plhs[0]) = (double) nrows;
                    return;
                }
            }
            else
                mexErrMsgTxt(mysql_error(conn));
        }
        ulong nrow = (ulong)mysql_num_rows(res), nfield = mysql_num_fields(res);
        //  If he didn't ask for any output (nlhs = 0),
        //     then display the output and return
        if (nlhs<1) {
            fancyprint(res);
            mysql_free_result(res);
            return;
        }
        //  If we are here, he wants output
        //  He must give exactly the right number of output arguments
        if (nlhs!=1)	{
            mysql_free_result(res);
            mexPrintf("You specified %d output arguments, and got %d columns of data\n", nlhs, nfield);
            mexErrMsgTxt("Must give one output argument");
        }
        //  Fix the column types to fix MySQL C API sloppiness
        MYSQL_FIELD*f = mysql_fetch_fields(res);
        //fix_types(f, res);
        //  Create field list
        const char **fieldnames=(const char **)mxMalloc(nfield*sizeof(char*));
        //
        for (ulong j=0; j<nfield;j++)
        {
            fieldnames[j]=f[j].name;
        }
        //  Create the Matlab structure for output
        
        if (!(plhs[0]=mxCreateStructMatrix(1,1,nfield,fieldnames)))
        {
            mysql_free_result(res);
            mexErrMsgTxt("Unable to create structure for output");
        }
        mxFree(fieldnames);
        // Create temporary mxArray pointer array
        mxArray* tmpPr;        
        //  Create the Matlab arrays for output        
        double**pr = (double**) mxMalloc(nfield*sizeof(double*));
        //    
        for (ulong j = 0; j<nfield; j++) {
            if (can_convert(f[j].type)) {
                if (!(tmpPr = mxCreateDoubleMatrix(nrow, 1, mxREAL))) {
                    mysql_free_result(res);
                    mexErrMsgTxt("Unable to create numeric matrix for output");
                }
                pr[j] = mxGetPr(tmpPr);
            }
            else if (f[j].type==FIELD_TYPE_BLOB) {
                if (!(tmpPr = mxCreateCellMatrix(nrow, 1))) {
                    mysql_free_result(res);
                    mexErrMsgTxt("Unable to create cell matrix for output");
                }
                pr[j]=NULL;
            }
            else {
                if (!(tmpPr = mxCreateCellMatrix(nrow, 1))) {
                    mysql_free_result(res);
                    mexErrMsgTxt("Unable to create cell matrix for output");
                }
                pr[j]=NULL;
            }
            mxSetField(plhs[0],0,f[j].name,tmpPr);
        }
        //  Load the data into the cells
        mysql_data_seek(res, 0);
        for (ulong i = 0; i<nrow; i++) {
            MYSQL_ROW row = mysql_fetch_row(res);
            if (!row) {
                mexPrintf("Scanning row %d for data extraction\n", i+1);
                mexErrMsgTxt("Internal error:  Failed to get a row");
            }
            ulong*p_lengths = mysql_fetch_lengths(res);
            for (ulong j = 0; j<nfield; j++) {
                if (can_convert(f[j].type)) {
                    pr[j][i] = field2num(row[j], f[j].type);
                }
                else if (f[j].type==FIELD_TYPE_BLOB) 
                {
                    tmpPr=mxGetField(plhs[0],0,f[j].name);
                    mxSetCell(tmpPr, i, deserialize(row[j], p_lengths[j]));
                }
                else {
                    tmpPr=mxGetField(plhs[0],0,f[j].name);
                    mxArray*c = mxCreateString(row[j]);
                    mxSetCell(tmpPr, i, c);
                }
            }
        }
        mxFree(pr);
        mysql_free_result(res);
    }
    else {
        mexPrintf("Unknown query type q = %d\n", q);
        mexErrMsgTxt("Internal code error");
    }
}
/*************************************SERIALIZE******************************************/
// Serialize a structure array
char* serializeStruct(ulong &rnBytes, const mxArray* rpArray, const char* rpArg, const bool rhead) {
    // number of dimensions
    int ndims = mxGetNumberOfDimensions(rpArray);
    rnBytes = 1+(1+ndims)*sizeof(int);
    if (rhead)
        rnBytes += (ulong)LEN_ID_MATLAB+1;
    // get number of fields and their names
    const int nfields = mxGetNumberOfFields(rpArray);
    rnBytes += sizeof(int);
    const char** pname = (const char**)mxCalloc(nfields, sizeof(char*));
    int* fnl = (int*)mxCalloc(nfields, sizeof(int));
    for (int i = 0; i<nfields; i++) {
        pname[i] = mxGetFieldNameByNumber(rpArray, i);
        fnl[i] = (int)strlen(pname[i]);
        rnBytes += fnl[i]+1;
    }
    // loop through each element of the array and serialize it
    const int nelts = mxGetNumberOfElements(rpArray);
    char** pser = (char**)mxCalloc(nelts*nfields, sizeof(char*));
    ulong* plen = (ulong*)mxCalloc(nelts*nfields, sizeof(ulong));
    const mxArray* pf;
    unsigned li = 0;
    for (int i = 0; i<nelts; i++)
        for (int j = 0; j<nfields; j++) {
            bool isWeCreatedAnEmptyMatrix=false;
			mxArray* fakeEmptyMatrix;
			pf = mxGetFieldByNumber(rpArray, i, j);
			if ((pf==NULL))
			{
				// sometimes, matlab returns null, for good and empty
				// elements. the reasons for this are unknown
				// so we create an empty double matrix, to overcome it
				fakeEmptyMatrix=mxCreateDoubleMatrix((mwSize)0,(mwSize)0,mxREAL);
				pf=fakeEmptyMatrix;
				isWeCreatedAnEmptyMatrix=true;
			}
			if (mxIsNumeric(pf)||mxIsLogical(pf)||mxIsChar(pf))
					pser[li] = serializeArray(plen[li], pf, rpArg, false);
				else if (mxIsStruct(pf))
					pser[li] = serializeStruct(plen[li], pf, rpArg, false);
				else if (mxIsCell(pf))
					pser[li] = serializeCell(plen[li], pf, rpArg, false);
				else
					mexErrMsgTxt("Unknow type");
			rnBytes += plen[li]+sizeof(int);
            li++;
			if (isWeCreatedAnEmptyMatrix)
			{
				mxDestroyArray(fakeEmptyMatrix);
			}
        }
    // put together the partial serializations
    char* p_serial = (char*)mxCalloc(rnBytes, sizeof(char));
    char* p_work = p_serial;
    if (rhead) {
        memcpy(p_work, ID_MATLAB, sizeof(char)*LEN_ID_MATLAB);
        p_work += LEN_ID_MATLAB;
        *(p_work++) = 0;
    }
    *(p_work++) = ID_STRUCT;
    // number of dimensions
    memcpy(p_work, (const char*)&ndims, sizeof(int));
    p_work += sizeof(int);
    // length of each dimension
    const int* pdims = mxGetDimensions(rpArray);
    memcpy(p_work, (const char*)pdims, sizeof(const int)*ndims);
    p_work += sizeof(const int)*ndims;
    // number of field
    memcpy(p_work, (const char*)&nfields, sizeof(const int));
    p_work += sizeof(const int);
    // field names
    for (int i = 0; i<nfields; i++) {
        memcpy(p_work, (const char*)pname[i], fnl[i]);
        p_work += fnl[i];
        *(p_work++) = 0;
    }
    // partial serialization
    int plenScalar;
    for (int i = 0; i<nfields*nelts; i++) {
        plenScalar=(int)plen[i];
        memcpy(p_work, (const char*)&plenScalar, sizeof(int));
        p_work  +=  sizeof(int);
        memcpy(p_work, (const char*)pser[i], plenScalar);
        p_work  +=  plen[i];
    }
    if (p_work-p_serial!=rnBytes)
        // should never arrives here
        mexErrMsgTxt("Oups, I haven't allocated enough memory, BEWARE");
    // free unused memory
    mxFree(pname);
    mxFree(fnl);
    for (int i = 1; i<nelts; i++)
        mxFree(pser[i]);
    mxFree(plen);
    return p_serial;
}
// Serialize a cell
char* serializeCell(ulong &rnBytes, const mxArray* rpArray, const char* rpArg, const bool rhead) {
    // number of dimensions
    const int ndims = mxGetNumberOfDimensions(rpArray);
    rnBytes = 1+(1+ndims)*sizeof(int);
    if (rhead)
        rnBytes += (u_long)LEN_ID_MATLAB+1;
    // loop through each element of the cell array and serialize it
    const int nelts = mxGetNumberOfElements(rpArray);
    char** pser = (char**)mxCalloc(nelts, sizeof(char*));
    ulong* plen = (ulong*)mxCalloc(nelts, sizeof(ulong));
    const mxArray* pf;
    for (int i = 0; i<nelts; i++){
        bool isWeCreatedAnEmptyMatrix=false;
		mxArray* fakeEmptyMatrix;
		pf = mxGetCell(rpArray, i);
		if (pf==NULL) 
		{	// sometimes pf null on empty elements,
			// we create a fake empty double element
			fakeEmptyMatrix=mxCreateDoubleMatrix((mwSize)0,(mwSize)0,mxREAL);
			pf=fakeEmptyMatrix;
			isWeCreatedAnEmptyMatrix=true;
		}
        if (mxIsNumeric(pf)||mxIsLogical(pf)||mxIsChar(pf))
            pser[i] = serializeArray(plen[i], pf, rpArg, false);
        else if (mxIsStruct(pf))
            pser[i] = serializeStruct(plen[i], pf, rpArg, false);
        else if (mxIsCell(pf))
            pser[i] = serializeCell(plen[i], pf, rpArg, false);
        else
            mexErrMsgTxt("unsupported type");
		rnBytes += plen[i]+sizeof(int);
		if (isWeCreatedAnEmptyMatrix)
		{
			mxDestroyArray(fakeEmptyMatrix);
		}
    }
    // put together the partial serializations
    char* p_serial = (char*)mxCalloc(rnBytes, sizeof(char));
    char* p_work = p_serial;
    if (rhead) {
        memcpy(p_work, ID_MATLAB, sizeof(char)*LEN_ID_MATLAB);
        p_work += LEN_ID_MATLAB;
        *(p_work++) = 0;
    }
    *(p_work++) = ID_CELL;
    // number of dimensions
    memcpy(p_work, (const char*)&ndims, sizeof(int));
    p_work += sizeof(int);
    // length of each dimension
    const int* pdims = mxGetDimensions(rpArray);
    memcpy(p_work, (const char*)pdims, sizeof(const int)*ndims);
    p_work += sizeof(const int)*ndims;
    // partial serialization
    int plenScalar;
    for (int i = 0; i<nelts; i++) {
        plenScalar=(int)plen[i];
        memcpy(p_work, (const char*)&plenScalar, sizeof(int));
        p_work  +=  sizeof(int);
        memcpy(p_work, (const char*)pser[i], plenScalar);
        p_work  +=  plen[i];
    }
    if (p_work-p_serial!=rnBytes)
        // should never arrives here
        mexErrMsgTxt("Oups, I haven't allocated enough memory, BEWARE");
    // free unused memory
    for (int i = 1; i<nelts; i++)
        mxFree(pser[i]);
    mxFree(plen);
    return p_serial;
}
// Serialize a numeric or a character array
char* serializeArray(ulong &rnBytes, const mxArray*rpArray, const char*rpArg, const bool rhead) {
    // compute number of bytes and allocate memory
    int ndims = mxGetNumberOfDimensions(rpArray);
    mxComplexity complexity = mxIsComplex(rpArray)? mxCOMPLEX : mxREAL;
    ulong n_bytes_data;
    if (complexity!=mxCOMPLEX)
        n_bytes_data = mxGetElementSize(rpArray)*mxGetNumberOfElements(rpArray);
    else
        n_bytes_data = mxGetElementSize(rpArray)*mxGetNumberOfElements(rpArray)*2;
    ulong nb = 1+sizeof(int)+ndims*sizeof(int)+sizeof(mxClassID)+sizeof(mxComplexity)+n_bytes_data;
    if (rhead)
        nb += (u_long)LEN_ID_MATLAB+1;
    char*p_serial = (char*)mxCalloc(nb, sizeof(char));
    char*p_work = p_serial;
    if (!p_serial)
        mexErrMsgTxt("Not enough memory");
    /***SERIALIZE***/
    // the first byte is set to 'A' for Array
    if (rhead) {
        memcpy(p_work, ID_MATLAB, sizeof(char)*LEN_ID_MATLAB);
        p_work += LEN_ID_MATLAB;
        *(p_work++) = 0;
    }
    *(p_work++) = ID_ARRAY;
    // number of dimensions
    memcpy(p_work, (const char*)&ndims, sizeof(int));
    p_work += sizeof(int);
    // length of each dimension
    const int*pdims = mxGetDimensions(rpArray);
    memcpy(p_work, (const char*)pdims, sizeof(const int)*ndims);
    p_work += sizeof(const int)*ndims;
    // classID
    mxClassID class_id = mxGetClassID(rpArray);
    memcpy(p_work, (const char*)&class_id, sizeof(mxClassID));
    p_work += sizeof(mxClassID);
    // complexity
    memcpy(p_work, (const char*)&complexity, sizeof(mxComplexity));
    p_work += sizeof(mxComplexity);
    // data
    if (complexity==mxCOMPLEX) {
        // real part data
        const char*pdata = (const char*)mxGetData(rpArray);
        memcpy(p_work, pdata, n_bytes_data/2);
        p_work += n_bytes_data/2;
        // imaginary part data (if exists)
        pdata = (const char*)mxGetImagData(rpArray);
        memcpy(p_work, pdata, n_bytes_data/2);
        p_work += n_bytes_data/2;
    }
    else{
        // real part data
        const char*pdata = (const char*)mxGetData(rpArray);
        memcpy(p_work, pdata, n_bytes_data);
        p_work += n_bytes_data;
    }
    rnBytes = (ulong)(p_work-p_serial);
    if (rnBytes > nb)
        // should never arrives here
        mexErrMsgTxt("Oups, I haven't allocated enough memory, BEWARE");
    return p_serial;
}
// Serialize a matlab string
char* serializeString(ulong &rnBytes, const mxArray*rpArray, const char*rpArg, const bool) {
    const int* pdims = mxGetDimensions(rpArray);
    const unsigned length = std::max<int>(pdims[0], pdims[1]);
    const unsigned n_dims = mxGetNumberOfDimensions(rpArray);
    char* p_buf = NULL;
    if (mxIsEmpty(rpArray)) {
        p_buf = (char*)mxCalloc(2, sizeof(char));
        p_buf[0] = ' ';
        p_buf[1] = 0;
        rnBytes = 1;
    }
    else if (mxIsChar(rpArray)) {
        if (n_dims!=2||!(pdims[0]==1||pdims[1]==1))
            mexErrMsgTxt("String placeholders only accept CHAR 1-by-M arrays or M-by-1!");
        // matlab string
        p_buf = (char*)mxCalloc(length+1, sizeof(char));
        mxGetString(rpArray, p_buf, length+1);
        rnBytes = length;
    }
    else if (mxIsNumeric(rpArray)||mxIsLogical(rpArray)) {
        // matlab scalar
        p_buf = (char*)mxCalloc(512, sizeof(char));
        if (length!=1)
            mexErrMsgTxt("Only scalars are accepted!");
        if (mxIsComplex(rpArray))
            mexErrMsgTxt("Complex scalars are not supported");
        double value;
        bool conv_to_int = true;
        if (mxIsDouble(rpArray)) {
            value = *((double*)mxGetData(rpArray));
            conv_to_int = false;
        }
        else if (mxIsSingle(rpArray)) {
            value = *((float*)mxGetData(rpArray));
            conv_to_int = false;
        }
        else if (mxIsInt8(rpArray))
            value = *((char*)mxGetData(rpArray));
        else if (mxIsInt32(rpArray))
            value = *((int*)mxGetData(rpArray));
        else if (mxIsUint8(rpArray))
            value = *((unsigned char*)mxGetData(rpArray));
        else if (mxIsUint32(rpArray))
            value = *((unsigned*)mxGetData(rpArray));
        else if (mxIsLogical(rpArray))
            value = *(char*)mxGetData(rpArray)!=0;
        else
            mexErrMsgTxt("Unsupported type, the supported types are (U)INT 8 and 32, FLOAT and DOUBLE");
        if (*rpArg=='i') {
            conv_to_int = true;
            value = floor(value);
        }
        if (conv_to_int)
            rnBytes = sprintf(p_buf, "%d", (int)value);
        else {
            int n = atoi(rpArg);
            if (n==0)
                rnBytes = sprintf(p_buf, "%.6f", value);
            else
                rnBytes = sprintf(p_buf, "%.*f", n, value);
        }
    }
    else
        mexErrMsgTxt("Only char string and non-complex scalar are supported");
    return p_buf;
}
// serialize binary data, given as a UINT8 vector
char* serializeBinary(ulong &rnBytes, const mxArray*rpArray, const char*rpArg, const bool) {
    rnBytes = mxGetNumberOfElements(rpArray);
    char*p_buf = (char*)mxCalloc(rnBytes, sizeof(char));
    memcpy(p_buf, (const char*)mxGetData(rpArray), rnBytes);
    return p_buf;
}
// Serialize a file, the filename is given
char* serializeFile(ulong &rnBytes, const mxArray*rpArray, const char*rpArg, const bool) {
    const int*pdims = mxGetDimensions(rpArray);
    const unsigned length = std::max<int>(pdims[0], pdims[1]);
    // filename provided
    char*p_name = (char*)mxCalloc(length+1, sizeof(char));
    mxGetString(rpArray, p_name, length+1);
    FILE*pf = fopen(p_name, "rb");
    mxFree(p_name);
    if (!pf)
        mexErrMsgTxt("Could not open file");
    rnBytes = file_length(pf);
    char*p_buf = (char*)mxCalloc(rnBytes, sizeof(char));
    fread(p_buf, sizeof(char), rnBytes, pf);
    fclose(pf);
    return p_buf;
}
/*************************************DESERIALIZE******************************************/
// deserialize a numeric or char Array
mxArray* deserializeArray(const char*rpSerial, const ulong rlength) {
    const char*p_serial = rpSerial;
    // number of dimensions
    int ndims;
    memcpy((char*)&ndims, p_serial, sizeof(int));
    p_serial += sizeof(int);
    // length of each dimension
    int*pdims = (int*)mxCalloc(ndims, sizeof(int));
    memcpy((char*)pdims, p_serial, sizeof(int)*ndims);
    p_serial += sizeof(int)*ndims;
    // classID
    mxClassID class_id;
    memcpy((char*)&class_id, p_serial, sizeof(mxClassID));
    p_serial += sizeof(mxClassID);
    // complexity
    mxComplexity complexity;
    memcpy((char*)&complexity, p_serial, sizeof(mxComplexity));
    p_serial += sizeof(mxComplexity);
    // create the numerice array
    mxArray*p_array = mxCreateNumericArray(ndims, pdims, class_id, complexity);
    mxFree(pdims);
    // real part data
    char*pdata = (char*)mxGetData(p_array);
    size_t nbytes_data_data = mxGetElementSize(p_array)*mxGetNumberOfElements(p_array);
    memcpy((char*)pdata, p_serial, nbytes_data_data);
    p_serial += nbytes_data_data;
    // imaginary part data (if exists)
    if (complexity==mxCOMPLEX) {
        pdata = (char*)mxGetImagData(p_array);
        memcpy((char*)pdata, p_serial, nbytes_data_data);
    }
    return p_array;
}
// deserialize struct Array
mxArray* deserializeStruct(const char* rpSerial, const ulong rlength){
    const char*p_serial = rpSerial;
    // number of dimensions
    int ndims;
    memcpy((char*)&ndims, p_serial, sizeof(int));
    p_serial += sizeof(int);
    // length of each dimension
    int* pdims = (int*)mxCalloc(ndims, sizeof(int));
    memcpy((char*)pdims, p_serial, sizeof(int)*ndims);
    p_serial += sizeof(int)*ndims;
    int nelts = 1;
    for (int i = 0; i<ndims; i++)
        nelts *= pdims[i];
    // number of field
    int nfields;
    memcpy(&nfields, (const char*)p_serial, sizeof(const int));
    p_serial += sizeof(const int);
    // field names
    const char** field_names = (const char**)mxCalloc(nfields, sizeof(char*));
    for (int i = 0; i<nfields; i++) {
        field_names[i] = p_serial;
        p_serial += strlen(p_serial)+1;
    }
    // create structure array
    mxArray* ps = mxCreateStructArray(ndims, pdims, nfields, field_names);
    // partial serialization
    for (int i = 0; i<nelts; i++)
        for (int j = 0; j<nfields; j++) {
            int nbytes1;
            ulong nbytes;
            memcpy(&nbytes1, p_serial, sizeof(int));
            nbytes=(long)nbytes1;
            p_serial  +=  sizeof(int);
            mxArray* pf;
            if (*p_serial==ID_ARRAY)
                pf = deserializeArray(p_serial+1, nbytes);
            else if (*p_serial==ID_STRUCT)
                pf = deserializeStruct(p_serial+1, nbytes);
            else if (*p_serial==ID_CELL)
                pf = deserializeCell(p_serial+1, nbytes);
            else
                mexErrMsgTxt("Unknow type");
            mxSetFieldByNumber(ps, i, j, pf);
            p_serial += nbytes;
        }
    // free unused memory
    mxFree(field_names);
    return ps;
}
// deserialize cell array
mxArray* deserializeCell(const char* rpSerial, const ulong rlength){
    const char* p_serial = rpSerial;
    // number of dimensions
    int ndims;
    memcpy((char*)&ndims, p_serial, sizeof(int));
    p_serial += sizeof(int);
    // length of each dimension
    int* pdims = (int*)mxCalloc(ndims, sizeof(int));
    memcpy((char*)pdims, p_serial, sizeof(int)*ndims);
    p_serial += sizeof(int)*ndims;
    int nelts = 1;
    for (int i = 0; i<ndims; i++)
        nelts *= pdims[i];
    // create structure array
    mxArray* ps = mxCreateCellArray(ndims, pdims);
    // partial serialization
    for (int i = 0; i<nelts; i++){
        ulong nbytes;
        int nbytes1;
        memcpy(&nbytes1, p_serial, sizeof(int));
        nbytes=(long)nbytes1;
        p_serial += sizeof(int);
        mxArray* pf;
        if (*p_serial==ID_ARRAY)
            pf = deserializeArray(p_serial+1, nbytes);
        else if (*p_serial==ID_STRUCT)
            pf = deserializeStruct(p_serial+1, nbytes);
        else if (*p_serial==ID_CELL)
            pf = deserializeCell(p_serial+1, nbytes);
        else
            mexErrMsgTxt("Unknow type");
        mxSetCell(ps, i, pf);
        p_serial += nbytes;
    }
    return ps;
}
/*
 *   returns the length of a file in bytes, and replace the file pointer to its state before calling
 */
int file_length(FILE*f) {
    int pos;
    int end;
    pos = ftell(f);
    fseek(f, 0, SEEK_END);
    end = ftell(f);
    fseek(f, pos, SEEK_SET);
    return end;
}
