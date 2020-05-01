/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2019 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: charmm_file.c,v $
 *      $Author: johns $        $Locker:  $             $State: Exp $
 *      $Revision: 1.4 $      $Date: 2019/07/24 03:36:56 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *  
 ***************************************************************************/

#include <stdio.h>
#include <ctype.h>
#include "charmm_file.h"

/*
   In addition to stream from which to read data, takes token buffer,
   max number of tokens, string buffer, and length of buffer.
   Strips comments and ignores lines containing no tokens.
   Reports legend lines (starting with *) with null first token.
   Returns number of tokens read or zero for end of file.
*/

int charmm_get_tokens(char **tok, int toklen,
			char *sbuf, int sbuflen,
			char *lbuf, int *lineno,
			FILE *stream, int all_caps) {

  int ntok;
  int fullline;
  char *s;
  char c2[2];

  ntok = 0;
  fullline = 0; /* Make compiler happy */
  while ( ! ntok ) {
    s = fgets(sbuf, sbuflen, stream);
    if ( ! s ) return 0;  /* EOF */
    if ( lbuf ) {
      char *l;
      for ( s = sbuf, l = lbuf; *s; ++s, ++l ) {
        *l = *s;
      }
      *l = *s;
    }
    if ( lineno ) ++(*lineno);
    for ( s = sbuf; *s; ++s ) {
      fullline = ( *s == '\n' );
    }
    if ( ! fullline ) do {
      s = fgets(c2, 2, stream);
    } while ( s && *s && *s != '\n' );
    s = sbuf;
    while ( *s && isspace(*s) ) ++s;
    if ( *s == '*' ) {
      *s = 0;
      tok[ntok] = s;  ++ntok;
      ++s;
      tok[ntok] = s;  ++ntok;
      while ( s && *s && *s != '\n' ) ++s;
      *s = 0;
      break;
    }
    if ( *s == '!' || *s == '\n' ) *s = 0;
    while ( *s ) {
      if ( ntok < toklen ) { tok[ntok] = s; ++ntok; }  /* in a token */
      while ( *s && *s != '!' && *s != '\n' && ! isspace(*s) ) {
        if ( all_caps ) *s = toupper(*s);
        ++s;
      }
      if ( *s == '!' || *s == '\n' ) *s = 0;
      if ( ! *s ) break;  /* no more tokens found on this line */
      *s = 0;  ++s;
      while ( *s && isspace(*s) ) ++s;
      if ( *s == '!' || *s == '\n' ) *s = 0;
    }
  }
  return ntok;

}

