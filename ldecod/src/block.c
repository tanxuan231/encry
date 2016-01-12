
/*!
 ***********************************************************************
 *  \file
 *      block.c
 *
 *  \brief
 *      Block functions
 *
 *  \author
 *      Main contributors (see contributors.h for copyright, address and affiliation details)
 *      - Inge Lille-Langoy          <inge.lille-langoy@telenor.com>
 *      - Rickard Sjoberg            <rickard.sjoberg@era.ericsson.se>
 ***********************************************************************
 */
#include "global.h"
#include "image.h"
#include "mb_access.h"
#include "quant.h"
#include "memalloc.h"

/*!
 ***********************************************************************
 * \brief
 *    Luma DC inverse transform
 ***********************************************************************
 */ 
//是否可以去掉
#if 1
void itrans_2(Macroblock *currMB,    //!< current macroblock
              ColorPlane pl)         //!< used color plane
{
  Slice *currSlice = currMB->p_Slice;
  VideoParameters *p_Vid = currMB->p_Vid;
  int j;

  int transform_pl = (p_Vid->separate_colour_plane_flag != 0) ? PLANE_Y : pl;
  int **cof = currSlice->cof[transform_pl];
  int qp_scaled = currMB->qp_scaled[transform_pl];

  int qp_per = p_Vid->qp_per_matrix[ qp_scaled ];
  int qp_rem = p_Vid->qp_rem_matrix[ qp_scaled ];      

  int invLevelScale = currSlice->InvLevelScale4x4_Intra[pl][qp_rem][0][0];
  int **M4;
  get_mem2Dint(&M4, BLOCK_SIZE, BLOCK_SIZE);
  
  // horizontal
  for (j=0; j < 4;++j) 
  {
    M4[j][0]=cof[j<<2][0];
    M4[j][1]=cof[j<<2][4];
    M4[j][2]=cof[j<<2][8];
    M4[j][3]=cof[j<<2][12];
  }

  ihadamard4x4(M4, M4);

  // vertical
  for (j=0; j < 4;++j) 
  {
    cof[j<<2][0]  = rshift_rnd((( M4[j][0] * invLevelScale) << qp_per), 6);
    cof[j<<2][4]  = rshift_rnd((( M4[j][1] * invLevelScale) << qp_per), 6);
    cof[j<<2][8]  = rshift_rnd((( M4[j][2] * invLevelScale) << qp_per), 6);
    cof[j<<2][12] = rshift_rnd((( M4[j][3] * invLevelScale) << qp_per), 6);
  }

  free_mem2Dint(M4);
}
#endif

void itrans_sp(Macroblock *currMB,   //!< current macroblock
               ColorPlane pl,        //!< used color plane
               int ioff,             //!< index to 4x4 block
               int joff)             //!< index to 4x4 block
{
  VideoParameters *p_Vid = currMB->p_Vid;
  Slice *currSlice = currMB->p_Slice;
  int i,j;  
  int ilev, icof;

  int qp = (currSlice->slice_type == SI_SLICE) ? currSlice->qs : currSlice->qp;
  int qp_per = p_Vid->qp_per_matrix[ qp ];
  int qp_rem = p_Vid->qp_rem_matrix[ qp ];

  int qp_per_sp = p_Vid->qp_per_matrix[ currSlice->qs ];
  int qp_rem_sp = p_Vid->qp_rem_matrix[ currSlice->qs ];
  int q_bits_sp = Q_BITS + qp_per_sp;

  imgpel **mb_pred = currSlice->mb_pred[pl];
  imgpel **mb_rec  = currSlice->mb_rec[pl];
  int    **mb_rres = currSlice->mb_rres[pl];
  int    **cof     = currSlice->cof[pl];
  int max_imgpel_value = p_Vid->max_pel_value_comp[pl];

  const int (*InvLevelScale4x4)  [4] = dequant_coef[qp_rem];
  const int (*InvLevelScale4x4SP)[4] = dequant_coef[qp_rem_sp];  
  int **PBlock;  

  get_mem2Dint(&PBlock, MB_BLOCK_SIZE, MB_BLOCK_SIZE);

  for (j=0; j< BLOCK_SIZE; ++j)
  {
    PBlock[j][0] = mb_pred[j+joff][ioff    ];
    PBlock[j][1] = mb_pred[j+joff][ioff + 1];
    PBlock[j][2] = mb_pred[j+joff][ioff + 2];
    PBlock[j][3] = mb_pred[j+joff][ioff + 3];
  }

  forward4x4(PBlock, PBlock, 0, 0);

  if(currSlice->sp_switch || currSlice->slice_type==SI_SLICE)
  {    
    for (j=0;j<BLOCK_SIZE;++j)
    {
      for (i=0;i<BLOCK_SIZE;++i)
      {
        // recovering coefficient since they are already dequantized earlier
        icof = (cof[joff + j][ioff + i] >> qp_per) / InvLevelScale4x4[j][i];
        //icof = ((cof[joff + j][ioff + i] * quant_coef[qp_rem][j][i])>> (qp_per + 15)) ;
        // icof  = rshift_rnd_sf(cof[joff + j][ioff + i] * quant_coef[qp_rem][j][i], qp_per + 15);
        ilev  = rshift_rnd_sf(iabs(PBlock[j][i]) * quant_coef[qp_rem_sp][j][i], q_bits_sp);
        ilev  = isignab(ilev, PBlock[j][i]) + icof;
        cof[joff + j][ioff + i] = ilev * InvLevelScale4x4SP[j][i] << qp_per_sp;
      }
    }
  }
  else
  {
    for (j=0;j<BLOCK_SIZE;++j)
    {
      for (i=0;i<BLOCK_SIZE;++i)
      {
        // recovering coefficient since they are already dequantized earlier
        icof = (cof[joff + j][ioff + i] >> qp_per) / InvLevelScale4x4[j][i];
        //icof = cof[joff + j][ioff + i];
        //icof  = rshift_rnd_sf(cof[joff + j][ioff + i] * quant_coef[qp_rem][j][i], qp_per + 15);
        ilev = PBlock[j][i] + ((icof * InvLevelScale4x4[j][i] * A[j][i] <<  qp_per) >> 6);
        ilev  = isign(ilev) * rshift_rnd_sf(iabs(ilev) * quant_coef[qp_rem_sp][j][i], q_bits_sp);
        //cof[joff + j][ioff + i] = ilev * InvLevelScale4x4SP[j][i] << qp_per_sp;
        cof[joff + j][ioff + i] = ilev * InvLevelScale4x4SP[j][i] << qp_per_sp;
      }
    }
  }

  inverse4x4(cof, mb_rres, joff, ioff);

  for (j=joff; j<joff +BLOCK_SIZE;++j)
  {
    mb_rec[j][ioff   ] = (imgpel) iClip1(max_imgpel_value,rshift_rnd_sf(mb_rres[j][ioff   ], DQ_BITS));
    mb_rec[j][ioff+ 1] = (imgpel) iClip1(max_imgpel_value,rshift_rnd_sf(mb_rres[j][ioff+ 1], DQ_BITS));
    mb_rec[j][ioff+ 2] = (imgpel) iClip1(max_imgpel_value,rshift_rnd_sf(mb_rres[j][ioff+ 2], DQ_BITS));
    mb_rec[j][ioff+ 3] = (imgpel) iClip1(max_imgpel_value,rshift_rnd_sf(mb_rres[j][ioff+ 3], DQ_BITS));
  }  

  free_mem2Dint(PBlock);
}

/*!
 *************************************************************************************
 * \brief
 *    Copy ImgPel Data from one structure to another (8x8)
 *************************************************************************************
 */
void copy_image_data_8x8(imgpel  **imgBuf1, imgpel  **imgBuf2, int off1, int off2)
{  
  int j;
  for(j = 0; j < BLOCK_SIZE_8x8; j+=4)
  {
    memcpy((*imgBuf1++ + off1), (*imgBuf2++ + off2), BLOCK_SIZE_8x8 * sizeof (imgpel));
    memcpy((*imgBuf1++ + off1), (*imgBuf2++ + off2), BLOCK_SIZE_8x8 * sizeof (imgpel));
    memcpy((*imgBuf1++ + off1), (*imgBuf2++ + off2), BLOCK_SIZE_8x8 * sizeof (imgpel));
    memcpy((*imgBuf1++ + off1), (*imgBuf2++ + off2), BLOCK_SIZE_8x8 * sizeof (imgpel));
  }
}
