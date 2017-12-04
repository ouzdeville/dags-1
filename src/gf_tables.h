/*
 * gf_tables.h
 *
 *  Created on: Dec 1, 2017
 *      Author: vader
 */

#ifndef SRC_GF_TABLES_H_
#define SRC_GF_TABLES_H_

#include <stdint.h>
#include "types_def.h"


gf_t gf_log_sf[32];
gf_t gf_antilog_sf[32];

gf gf_log[1024];
gf gf_antilog[1024];

#endif /* SRC_GF_TABLES_H_ */
