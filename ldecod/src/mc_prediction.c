
/*!
 *************************************************************************************
 * \file mc_prediction.c
 *
 * \brief
 *    Functions for motion compensated prediction
 *
 * \author
 *      Main contributors (see contributors.h for copyright, 
 *                         address and affiliation details)
 *      - Alexis Michael Tourapis  <alexismt@ieee.org>
 *      - Chris Vogt
 *
 *************************************************************************************
 */
#include "global.h"
#include "mbuffer.h"
#include "mb_access.h"
#include "macroblock.h"
#include "memalloc.h"
#include "dec_statistics.h"

int allocate_pred_mem(Slice *currSlice)
{
  int alloc_size = 0;
  alloc_size += get_mem2Dpel(&currSlice->tmp_block_l0, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
  alloc_size += get_mem2Dpel(&currSlice->tmp_block_l1, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
  alloc_size += get_mem2Dpel(&currSlice->tmp_block_l2, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
  alloc_size += get_mem2Dpel(&currSlice->tmp_block_l3, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
  alloc_size += get_mem2Dint(&currSlice->tmp_res, MB_BLOCK_SIZE + 5, MB_BLOCK_SIZE + 5);
  return (alloc_size);
}

void free_pred_mem(Slice *currSlice)
{
  free_mem2Dint(currSlice->tmp_res);
  free_mem2Dpel(currSlice->tmp_block_l0);
  free_mem2Dpel(currSlice->tmp_block_l1);
  free_mem2Dpel(currSlice->tmp_block_l2);
  free_mem2Dpel(currSlice->tmp_block_l3);
}
