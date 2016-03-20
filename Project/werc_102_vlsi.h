/*************************************************
* File name   : werc_102_vlsi.h
* Description : Refine the assert method and data type of all the project
*               This head file should be contained by all the project
* Author      : WJF
* Date        :
* Compiler    : gcc-4.6.3
* Target      : Ubuntu
* History     :
*   <author>  <time>   <desc>
*     WJF    20130808  CREATE
*************************************************/


#ifndef _WERC_102_VLSI_H_
#define _WERC_102_VLSI_H_


#ifndef NDEBUG
#include    <assert.h>
#define H102_ASSERT(_exp)    assert(_exp)
#else
#define H102_ASSERT(_exp)   (void 0)
#endif


/** MACRO OF BYTE ORDER **/
#define SWAP16(x)   ((((x) <<  8) & 0xFF00) | (((x) >>  8) & 0xFF))
#define SWAP32(x)   ((((x) << 24) & 0xFF000000)  | \
	(((x) << 8) & 0xFF0000) | \
	(((x) >> 8) & 0xFF00) | \
	(((x) >> 24) & 0xFF))

#if INSTALL_ON_BIG_ENDIAN
#define ECEN_HTONS(x)     (x)
#define ECEN_NTOHS(x)     (x)
#define ECEN_HTONL(x)     (x)
#define ECEN_NTOHL(x)     (x)
#else /*_LITTLE_ENDIAN*/
#define ECEN_HTONS(x)     SWAP16(x)
#define ECEN_NTOHS(x)     SWAP16(x)
#define ECEN_HTONL(x)     SWAP32(x)
#define ECEN_NTOHL(x)     SWAP32(x)
#endif /*end of #ifdef _BIG_ENDIAN*/

/** MACRO OF BIT SET **/
#ifndef BIT_SET
#define BIT_SET(F, B)   ((F) |= (B))
#endif // BIT_SET

#ifndef BIT_RESET
#define BIT_RESET(F, B) ((F) &= ~(B))
#endif // BIT_RESET

#ifndef BIT_TEST
#define BIT_TEST(F, B)  ((F) & (B))
#endif // BIT_TEST

/** Basic Data Structure **/
typedef     void    VOID;
typedef     char    CHAR;
typedef     unsigned char   UINT8D;
typedef     unsigned short  UINT16D;
typedef     unsigned int    UINT32D;
typedef     unsigned long   UINT64D;
typedef     char    INT8D;
typedef     short   INT16D;
typedef     int     INT32D;
typedef     long    INT64D;
typedef     float   FLT32;
typedef     double  FLT64;

typedef     unsigned char   BOOLEAN;

#define     NO      0
#define     YES     1
#define     FALSE   0
#define     TRUE    1

#ifndef NULL
#define NULL ((VOID*)0)
#endif

#endif

