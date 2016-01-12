
/*!
 ***********************************************************************
 * \file image.c
 *
 * \brief
 *    Decode a Slice
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Inge Lille-Langoy               <inge.lille-langoy@telenor.com>
 *    - Rickard Sjoberg                 <rickard.sjoberg@era.ericsson.se>
 *    - Jani Lainema                    <jani.lainema@nokia.com>
 *    - Sebastian Purreiter             <sebastian.purreiter@mch.siemens.de>
 *    - Byeong-Moon Jeon                <jeonbm@lge.com>
 *    - Thomas Wedi                     <wedi@tnt.uni-hannover.de>
 *    - Gabi Blaettermann
 *    - Ye-Kui Wang                     <wyk@ieee.org>
 *    - Antti Hallapuro                 <antti.hallapuro@nokia.com>
 *    - Alexis Tourapis                 <alexismt@ieee.org>
 *    - Jill Boyce                      <jill.boyce@thomson.net>
 *    - Saurav K Bandyopadhyay          <saurav@ieee.org>
 *    - Zhenyu Wu                       <Zhenyu.Wu@thomson.net
 *    - Purvin Pandit                   <Purvin.Pandit@thomson.net>
 *
 ***********************************************************************
 */
#include <math.h>
#include <limits.h>

#include "global.h"
#include "image.h"
#include "fmo.h"
#include "annexb.h"
#include "nalu.h"
#include "parset.h"
#include "header.h"

#include "sei.h"
#include "mb_access.h"
#include "memalloc.h"
#include "macroblock.h"

#include "biaridecod.h"
#include "context_ini.h"
#include "cabac.h"
#include "vlc.h"
#include "quant.h"

#include "fast_memory.h"

extern int testEndian(void);

static inline void reset_mbs(Macroblock *currMB)
{
  currMB->slice_nr = -1; 
  currMB->ei_flag  =  1;
  currMB->dpl_flag =  0;
}

static void setup_buffers(VideoParameters *p_Vid, int layer_id)
{
  CodingParameters *cps = p_Vid->p_EncodePar[layer_id];
  int i;

  if(p_Vid->last_dec_layer_id != layer_id)
  {
    //p_Vid->imgY_ref = cps->imgY_ref;
    //p_Vid->imgUV_ref = cps->imgUV_ref;
    if(cps->separate_colour_plane_flag)
    {
     for( i=0; i<MAX_PLANE; i++ )
     {
       p_Vid->mb_data_JV[i] = cps->mb_data_JV[i];
       p_Vid->intra_block_JV[i] = cps->intra_block_JV[i];
       p_Vid->ipredmode_JV[i] = cps->ipredmode_JV[i];
       p_Vid->siblock_JV[i] = cps->siblock_JV[i];
     }
     p_Vid->mb_data = NULL;
     p_Vid->intra_block = NULL;
     p_Vid->ipredmode = NULL;
     p_Vid->siblock = NULL;
    }
    else
    {
      p_Vid->mb_data = cps->mb_data;
      p_Vid->intra_block = cps->intra_block;
      p_Vid->ipredmode = cps->ipredmode;
      p_Vid->siblock = cps->siblock;
    }
    p_Vid->PicPos = cps->PicPos;
    p_Vid->nz_coeff = cps->nz_coeff;
    p_Vid->qp_per_matrix = cps->qp_per_matrix;
    p_Vid->qp_rem_matrix = cps->qp_rem_matrix;
    p_Vid->oldFrameSizeInMbs = cps->oldFrameSizeInMbs;
    p_Vid->img2buf = cps->img2buf;
    p_Vid->last_dec_layer_id = layer_id;
  }
}

#if MVC_EXTENSION_ENABLE
static void init_mvc_picture(Slice *currSlice)
{
  int i;
  VideoParameters *p_Vid = currSlice->p_Vid;
  //DecodedPictureBuffer *p_Dpb = p_Vid->p_Dpb_layer[0];

  StorablePicture *p_pic = NULL;

#if 0
  // find BL reconstructed picture
  if (currSlice->structure  == FRAME)
  {
    for (i = 0; i < (int)p_Dpb->used_size/*size*/; i++)
    {
      FrameStore *fs = p_Dpb->fs[i];
      if ((fs->frame->view_id == 0) && (fs->frame->frame_poc == currSlice->framepoc))
      {
        p_pic = fs->frame;
        break;
      }
    }
  }
  else if (currSlice->structure  == TOP_FIELD)
  {
    for (i = 0; i < (int)p_Dpb->used_size/*size*/; i++)
    {
      FrameStore *fs = p_Dpb->fs[i];
      if ((fs->top_field->view_id == 0) && (fs->top_field->top_poc == currSlice->toppoc))
      {
        p_pic = fs->top_field;
        break;
      }
    }
  }
  else
  {
    for (i = 0; i < (int)p_Dpb->used_size/*size*/; i++)
    {
      FrameStore *fs = p_Dpb->fs[i];
      if ((fs->bottom_field->view_id == 0) && (fs->bottom_field->bottom_poc == currSlice->bottompoc))
      {
        p_pic = fs->bottom_field;
        break;
      }
    }
  }
#endif  
  if(!p_pic)
  {
    p_Vid->bFrameInit = 0;
  }
  else
  {
    //process_picture_in_dpb_s(p_Vid, p_pic);
    //store_proc_picture_in_dpb (currSlice->p_Dpb, clone_storable_picture(p_Vid, p_pic));
  }
}
#endif

/*!
 ************************************************************************
 * \brief
 *    Initializes the parameters for a new picture
 ************************************************************************
 */
//当读取到新一帧的第一个片时调用  (需要处理完p_Vid->dec_picture)
#if 1
static void init_picture(VideoParameters *p_Vid, Slice *currSlice, InputParameters *p_Inp)
{
  int i;
  int nplane;
  StorablePicture *dec_picture = NULL;
  seq_parameter_set_rbsp_t *active_sps = p_Vid->active_sps;
  //DecodedPictureBuffer *p_Dpb = currSlice->p_Dpb;

  p_Vid->PicHeightInMbs = p_Vid->FrameHeightInMbs / ( 1 + currSlice->field_pic_flag );
  p_Vid->PicSizeInMbs   = p_Vid->PicWidthInMbs * p_Vid->PicHeightInMbs;
  p_Vid->FrameSizeInMbs = p_Vid->PicWidthInMbs * p_Vid->FrameHeightInMbs;

  p_Vid->bFrameInit = 1;
  if (p_Vid->dec_picture) // && p_Vid->num_dec_mb == p_Vid->PicSizeInMbs)
  {
    // this may only happen on slice loss
    exit_picture(p_Vid, &p_Vid->dec_picture);
  }
  p_Vid->dpb_layer_id = currSlice->layer_id;
  //set buffers;
  setup_buffers(p_Vid, currSlice->layer_id);

  if (p_Vid->recovery_point)
    p_Vid->recovery_frame_num = (currSlice->frame_num + p_Vid->recovery_frame_cnt) % p_Vid->max_frame_num;

  if (currSlice->idr_flag)
    p_Vid->recovery_frame_num = currSlice->frame_num;

  if(currSlice->nal_reference_idc)
  {
    p_Vid->pre_frame_num = currSlice->frame_num;
  }

  //p_Vid->num_dec_mb = 0;

  //calculate POC
  //decode_poc(p_Vid, currSlice);

  //if (p_Vid->recovery_frame_num == (int) currSlice->frame_num && p_Vid->recovery_poc == 0x7fffffff)
    //p_Vid->recovery_poc = currSlice->framepoc;

  //if(currSlice->nal_reference_idc)
    //p_Vid->last_ref_pic_poc = currSlice->framepoc;

  //  dumppoc (p_Vid);

  if (currSlice->structure==FRAME ||currSlice->structure==TOP_FIELD)
  {
    gettime (&(p_Vid->start_time));             // start time
  }

  dec_picture = p_Vid->dec_picture = alloc_storable_picture (p_Vid, currSlice->structure, p_Vid->width, p_Vid->height, p_Vid->width_cr, p_Vid->height_cr, 1);
  //dec_picture->top_poc=currSlice->toppoc;
  //dec_picture->bottom_poc=currSlice->bottompoc;
  //dec_picture->frame_poc=currSlice->framepoc;
  dec_picture->qp = currSlice->qp;
  dec_picture->slice_qp_delta = currSlice->slice_qp_delta;
  dec_picture->chroma_qp_offset[0] = p_Vid->active_pps->chroma_qp_index_offset;
  dec_picture->chroma_qp_offset[1] = p_Vid->active_pps->second_chroma_qp_index_offset;
  dec_picture->iCodingType = currSlice->structure==FRAME? (currSlice->mb_aff_frame_flag? FRAME_MB_PAIR_CODING:FRAME_CODING): FIELD_CODING; //currSlice->slice_type;
  dec_picture->layer_id = currSlice->layer_id;
#if (MVC_EXTENSION_ENABLE)
  dec_picture->view_id         = currSlice->view_id;
  dec_picture->inter_view_flag = currSlice->inter_view_flag;
  dec_picture->anchor_pic_flag = currSlice->anchor_pic_flag;
  if (dec_picture->view_id == 1)
  {
    if((p_Vid->profile_idc == MVC_HIGH) || (p_Vid->profile_idc == STEREO_HIGH))
      init_mvc_picture(currSlice);
  }
#endif

  // reset all variables of the error concealment instance before decoding of every frame.
  // here the third parameter should, if perfectly, be equal to the number of slices per frame.
  // using little value is ok, the code will allocate more memory if the slice number is larger
#if (DISABLE_ERC == 0)
  //ercReset(p_Vid->erc_errorVar, p_Vid->PicSizeInMbs, p_Vid->PicSizeInMbs, dec_picture->size_x);
#endif
  //p_Vid->erc_mvperMB = 0;

  switch (currSlice->structure )
  {
  case TOP_FIELD:
    {
      //dec_picture->poc = currSlice->toppoc;
      p_Vid->number *= 2;
      break;
    }
  case BOTTOM_FIELD:
    {
      //dec_picture->poc = currSlice->bottompoc;
      p_Vid->number = p_Vid->number * 2 + 1;
      break;
    }
  case FRAME:
    {
      //dec_picture->poc = currSlice->framepoc;
      break;
    }
  default:
    error("p_Vid->structure not initialized", 235);
  }

  //p_Vid->current_slice_nr=0;

  if (p_Vid->type > SI_SLICE)
  {
    //set_ec_flag(p_Vid, SE_PTYPE);
    p_Vid->type = P_SLICE;  // concealed element
  }

  // CAVLC init
  if (p_Vid->active_pps->entropy_coding_mode_flag == (Boolean) CAVLC)
  {
    memset(p_Vid->nz_coeff[0][0][0], -1, p_Vid->PicSizeInMbs * 48 *sizeof(byte)); // 3 * 4 * 4
  }

  // Set the slice_nr member of each MB to -1, to ensure correct when packet loss occurs
  // TO set Macroblock Map (mark all MBs as 'have to be concealed')
  if( (p_Vid->separate_colour_plane_flag != 0) )
  {
    for( nplane=0; nplane<MAX_PLANE; ++nplane )
    {      
      Macroblock *currMB = p_Vid->mb_data_JV[nplane];
      char *intra_block = p_Vid->intra_block_JV[nplane];
      for(i=0; i<(int)p_Vid->PicSizeInMbs; ++i)
      {
        reset_mbs(currMB++);
      }
      fast_memset(p_Vid->ipredmode_JV[nplane][0], DC_PRED, 16 * p_Vid->FrameHeightInMbs * p_Vid->PicWidthInMbs * sizeof(char));
      if(p_Vid->active_pps->constrained_intra_pred_flag)
      {
        for (i=0; i<(int)p_Vid->PicSizeInMbs; ++i)
        {
          intra_block[i] = 1;
        }
      }
    }
  }
  else
  {
#if 0 //defined(OPENMP)
#pragma omp parallel for
    for(i=0; i<(int)p_Vid->PicSizeInMbs; ++i)
      reset_mbs(&p_Vid->mb_data[i]);
#else
    Macroblock *currMB = p_Vid->mb_data;
    for(i=0; i<(int)p_Vid->PicSizeInMbs; ++i)
      reset_mbs(currMB++);
#endif
    if(p_Vid->active_pps->constrained_intra_pred_flag)
    {
      for (i=0; i<(int)p_Vid->PicSizeInMbs; ++i)
      {
        p_Vid->intra_block[i] = 1;
      }
    }
    fast_memset(p_Vid->ipredmode[0], DC_PRED, 16 * p_Vid->FrameHeightInMbs * p_Vid->PicWidthInMbs * sizeof(char));
  }  

  dec_picture->slice_type = p_Vid->type;
  dec_picture->used_for_reference = (currSlice->nal_reference_idc != 0);
  dec_picture->idr_flag = currSlice->idr_flag;
  dec_picture->no_output_of_prior_pics_flag = currSlice->no_output_of_prior_pics_flag;
  dec_picture->long_term_reference_flag     = currSlice->long_term_reference_flag;
  dec_picture->adaptive_ref_pic_buffering_flag = currSlice->adaptive_ref_pic_buffering_flag;

  dec_picture->dec_ref_pic_marking_buffer = currSlice->dec_ref_pic_marking_buffer;
  currSlice->dec_ref_pic_marking_buffer   = NULL;

  dec_picture->mb_aff_frame_flag = currSlice->mb_aff_frame_flag;
  dec_picture->PicWidthInMbs     = p_Vid->PicWidthInMbs;

  p_Vid->get_mb_block_pos = dec_picture->mb_aff_frame_flag ? get_mb_block_pos_mbaff : get_mb_block_pos_normal;
  p_Vid->getNeighbour     = dec_picture->mb_aff_frame_flag ? getAffNeighbour : getNonAffNeighbour;

  dec_picture->pic_num   = currSlice->frame_num;
  dec_picture->frame_num = currSlice->frame_num;

  dec_picture->recovery_frame = (unsigned int) ((int) currSlice->frame_num == p_Vid->recovery_frame_num);

  dec_picture->coded_frame = (currSlice->structure==FRAME);

  dec_picture->chroma_format_idc = active_sps->chroma_format_idc;

  dec_picture->frame_mbs_only_flag = active_sps->frame_mbs_only_flag;
  dec_picture->frame_cropping_flag = active_sps->frame_cropping_flag;

  if (dec_picture->frame_cropping_flag)
  {
    dec_picture->frame_crop_left_offset   = active_sps->frame_crop_left_offset;
    dec_picture->frame_crop_right_offset  = active_sps->frame_crop_right_offset;
    dec_picture->frame_crop_top_offset    = active_sps->frame_crop_top_offset;
    dec_picture->frame_crop_bottom_offset = active_sps->frame_crop_bottom_offset;
  }

#if (ENABLE_OUTPUT_TONEMAPPING)
  // store the necessary tone mapping sei into StorablePicture structure
  if (p_Vid->seiToneMapping->seiHasTone_mapping)
  {
    int coded_data_bit_max = (1 << p_Vid->seiToneMapping->coded_data_bit_depth);
    dec_picture->seiHasTone_mapping    = 1;
    dec_picture->tone_mapping_model_id = p_Vid->seiToneMapping->model_id;
    dec_picture->tonemapped_bit_depth  = p_Vid->seiToneMapping->sei_bit_depth;
    dec_picture->tone_mapping_lut      = malloc(coded_data_bit_max * sizeof(int));
    if (NULL == dec_picture->tone_mapping_lut)
    {
      no_mem_exit("init_picture: tone_mapping_lut");
    }
    memcpy(dec_picture->tone_mapping_lut, p_Vid->seiToneMapping->lut, sizeof(imgpel) * coded_data_bit_max);
    update_tone_mapping_sei(p_Vid->seiToneMapping);
  }
  else
    dec_picture->seiHasTone_mapping = 0;
#endif

  if( (p_Vid->separate_colour_plane_flag != 0) )
  {
    //p_Vid->dec_picture_JV[0] = p_Vid->dec_picture;
    //p_Vid->dec_picture_JV[1] = alloc_storable_picture (p_Vid, (PictureStructure) currSlice->structure, p_Vid->width, p_Vid->height, p_Vid->width_cr, p_Vid->height_cr, 1);
    //copy_dec_picture_JV( p_Vid, p_Vid->dec_picture_JV[1], p_Vid->dec_picture_JV[0] );
    //p_Vid->dec_picture_JV[2] = alloc_storable_picture (p_Vid, (PictureStructure) currSlice->structure, p_Vid->width, p_Vid->height, p_Vid->width_cr, p_Vid->height_cr, 1);
    //copy_dec_picture_JV( p_Vid, p_Vid->dec_picture_JV[2], p_Vid->dec_picture_JV[0] );
  }
}
#endif

static void update_mbaff_macroblock_data(imgpel **cur_img, imgpel (*temp)[16], int x0, int width, int height)
{
  imgpel (*temp_evn)[16] = temp;
  imgpel (*temp_odd)[16] = temp + height; 
  imgpel **temp_img = cur_img;
  int y;

  for (y = 0; y < 2 * height; ++y)
    memcpy(*temp++, (*temp_img++ + x0), width * sizeof(imgpel));

  for (y = 0; y < height; ++y)
  {
    memcpy((*cur_img++ + x0), *temp_evn++, width * sizeof(imgpel));
    memcpy((*cur_img++ + x0), *temp_odd++, width * sizeof(imgpel));
  }
}

static void MbAffPostProc(VideoParameters *p_Vid)
{

  imgpel temp_buffer[32][16];

  StorablePicture *dec_picture = p_Vid->dec_picture;
  imgpel ** imgY  = dec_picture->imgY;
  imgpel ***imgUV = dec_picture->imgUV;

  short i, x0, y0;

  for (i=0; i<(int)dec_picture->PicSizeInMbs; i+=2)
  {
    if (dec_picture->motion.mb_field[i])
    {
      get_mb_pos(p_Vid, i, p_Vid->mb_size[IS_LUMA], &x0, &y0);
      update_mbaff_macroblock_data(imgY + y0, temp_buffer, x0, MB_BLOCK_SIZE, MB_BLOCK_SIZE);

      if (dec_picture->chroma_format_idc != YUV400)
      {
        x0 = (short) ((x0 * p_Vid->mb_cr_size_x) >> 4);
        y0 = (short) ((y0 * p_Vid->mb_cr_size_y) >> 4);

        update_mbaff_macroblock_data(imgUV[0] + y0, temp_buffer, x0, p_Vid->mb_cr_size_x, p_Vid->mb_cr_size_y);
        update_mbaff_macroblock_data(imgUV[1] + y0, temp_buffer, x0, p_Vid->mb_cr_size_x, p_Vid->mb_cr_size_y);
      }
    }
  }
}


static void init_picture_decoding(VideoParameters *p_Vid)
{
  Slice *pSlice = p_Vid->ppSliceList[0];
  int j, i, iDeblockMode=1;

  if(p_Vid->iSliceNumOfCurrPic >= MAX_NUM_SLICES)
  {
    error ("Maximum number of supported slices exceeded. \nPlease recompile with increased value for MAX_NUM_SLICES", 200);
  }

  if(p_Vid->pNextPPS->Valid && (int) p_Vid->pNextPPS->pic_parameter_set_id == pSlice->pic_parameter_set_id)
  {
    pic_parameter_set_rbsp_t tmpPPS;
    memcpy(&tmpPPS, &(p_Vid->PicParSet[pSlice->pic_parameter_set_id]), sizeof (pic_parameter_set_rbsp_t));
    (p_Vid->PicParSet[pSlice->pic_parameter_set_id]).slice_group_id = NULL;
    MakePPSavailable (p_Vid, p_Vid->pNextPPS->pic_parameter_set_id, p_Vid->pNextPPS);
    memcpy(p_Vid->pNextPPS, &tmpPPS, sizeof (pic_parameter_set_rbsp_t));
    tmpPPS.slice_group_id = NULL;
  }

  UseParameterSet (pSlice);
  if(pSlice->idr_flag)
    p_Vid->number=0;

  p_Vid->PicHeightInMbs = p_Vid->FrameHeightInMbs / ( 1 + pSlice->field_pic_flag );
  p_Vid->PicSizeInMbs   = p_Vid->PicWidthInMbs * p_Vid->PicHeightInMbs;
  p_Vid->FrameSizeInMbs = p_Vid->PicWidthInMbs * p_Vid->FrameHeightInMbs;
  p_Vid->structure = pSlice->structure;

  fmo_init (p_Vid, pSlice);

#if (MVC_EXTENSION_ENABLE)
  if((pSlice->layer_id>0) && (pSlice->svc_extension_flag == 0 && pSlice->NaluHeaderMVCExt.non_idr_flag == 0))
  {
   //idr_memory_management(p_Vid->p_Dpb_layer[pSlice->layer_id], p_Vid->dec_picture);
  }
  //update_ref_list(p_Vid->p_Dpb_layer[pSlice->view_id]);
  //update_ltref_list(p_Vid->p_Dpb_layer[pSlice->view_id]);
  //update_pic_num(pSlice);
  i = pSlice->view_id;
#else
  //update_pic_num(pSlice);
  i = 0;
#endif

  //init_Deblock(p_Vid, pSlice->mb_aff_frame_flag);	//需要去掉
#if 0
  //init mb_data;
  for(j=0; j<p_Vid->iSliceNumOfCurrPic; j++)
  {
    if(p_Vid->ppSliceList[j]->DFDisableIdc != 1)
      iDeblockMode=0;
#if (MVC_EXTENSION_ENABLE)
    assert(p_Vid->ppSliceList[j]->view_id == i);
#endif
  }
#endif  
}

void init_slice(VideoParameters *p_Vid, Slice *currSlice)
{
  int i;

  p_Vid->active_sps = currSlice->active_sps;
  p_Vid->active_pps = currSlice->active_pps;

  //currSlice->init_lists (currSlice);	//need fix

#if 0
#if (MVC_EXTENSION_ENABLE)
  if (currSlice->svc_extension_flag == 0 || currSlice->svc_extension_flag == 1)
    //reorder_lists_mvc (currSlice, currSlice->ThisPOC);
  else
    //reorder_lists (currSlice);

  if (currSlice->fs_listinterview0)
  {
    free(currSlice->fs_listinterview0);
    currSlice->fs_listinterview0 = NULL;
  }
  if (currSlice->fs_listinterview1)
  {
    free(currSlice->fs_listinterview1);
    currSlice->fs_listinterview1 = NULL;
  }
#else
  //reorder_lists (currSlice);
#endif
#endif

  if (currSlice->structure==FRAME)
  {
    //init_mbaff_lists(p_Vid, currSlice);
  }
  //p_Vid->recovery_point = 0;

#if 0
  // update reference flags and set current p_Vid->ref_flag
  if(!(currSlice->redundant_pic_cnt != 0 && p_Vid->previous_frame_num == currSlice->frame_num))
  {
    for(i=16;i>0;i--)
    {
      currSlice->ref_flag[i] = currSlice->ref_flag[i-1];
    }
  }
  currSlice->ref_flag[0] = currSlice->redundant_pic_cnt==0 ? p_Vid->Is_primary_correct : p_Vid->Is_redundant_correct;
  //p_Vid->previous_frame_num = currSlice->frame_num; //p_Vid->frame_num;
#endif

  if((currSlice->active_sps->chroma_format_idc==0)||(currSlice->active_sps->chroma_format_idc==3))
  {
    currSlice->linfo_cbp_intra = linfo_cbp_intra_other;
    currSlice->linfo_cbp_inter = linfo_cbp_inter_other;
  }
  else
  {
    currSlice->linfo_cbp_intra = linfo_cbp_intra_normal;
    currSlice->linfo_cbp_inter = linfo_cbp_inter_normal;
  }
}

void decode_slice(Slice *currSlice, int current_header)
{
  if (currSlice->active_pps->entropy_coding_mode_flag)
  {
    init_contexts  (currSlice);
    cabac_new_slice(currSlice);
  }

#if 0	//what's this ?
  if ( (currSlice->active_pps->weighted_bipred_idc > 0  && (currSlice->slice_type == B_SLICE)) || (currSlice->active_pps->weighted_pred_flag && currSlice->slice_type !=I_SLICE))
    fill_wp_params(currSlice);	//what's this? 
#endif 

  // decode main slice information
  if ((current_header == SOP || current_header == SOS) && currSlice->ei_flag == 0)
    decode_one_slice(currSlice);

  // setMB-Nr in case this slice was lost
  // if(currSlice->ei_flag)
  //   p_Vid->current_mb_nr = currSlice->last_mb_nr + 1;

}

#if 0
/*!
 ************************************************************************
 * \brief
 *    Error tracking: if current frame is lost or any reference frame of
 *                    current frame is lost, current frame is incorrect.
 ************************************************************************
 */
static void Error_tracking(VideoParameters *p_Vid, Slice *currSlice)
{
  int i;

  if(currSlice->redundant_pic_cnt == 0)
  {
    p_Vid->Is_primary_correct = p_Vid->Is_redundant_correct = 1;
  }

  if(currSlice->redundant_pic_cnt == 0 && p_Vid->type != I_SLICE)
  {
    for(i=0;i<currSlice->num_ref_idx_active[LIST_0];++i)
    {
      if(currSlice->ref_flag[i] == 0)  // any reference of primary slice is incorrect
      {
        p_Vid->Is_primary_correct = 0; // primary slice is incorrect
      }
    }
  }
  else if(currSlice->redundant_pic_cnt != 0 && p_Vid->type != I_SLICE)
  {
    if(currSlice->ref_flag[currSlice->redundant_slice_ref_idx] == 0)  // reference of redundant slice is incorrect
    {
      p_Vid->Is_redundant_correct = 0;  // redundant slice is incorrect
    }
  }
}
#endif

static void CopyPOC(Slice *pSlice0, Slice *currSlice)
{
  //currSlice->framepoc  = pSlice0->framepoc;
  //currSlice->toppoc    = pSlice0->toppoc;
  //currSlice->bottompoc = pSlice0->bottompoc;  
  //currSlice->ThisPOC   = pSlice0->ThisPOC;
}



/*!
 ***********************************************************************
 * \brief
 *    decodes one I- or P-frame
 *
 ***********************************************************************
 */
int decode_one_frame(DecoderParams *pDecoder)
{
  VideoParameters *p_Vid = pDecoder->p_Vid;
  InputParameters *p_Inp = p_Vid->p_Inp;
  int current_header, iRet;
  Slice *currSlice; // = p_Vid->currentSlice;
  Slice **ppSliceList = p_Vid->ppSliceList;
  int iSliceNo;
  
  //read one picture first;
  p_Vid->iSliceNumOfCurrPic=0;  //当前图像的片序号
  current_header=0;
  p_Vid->iNumOfSlicesDecoded=0;  //解码片的序号
  p_Vid->num_dec_mb = 0;
	
  if(p_Vid->newframe)
  {
    if(p_Vid->pNextPPS->Valid) //下一个PPS有效
    {
      //assert((int) p_Vid->pNextPPS->pic_parameter_set_id == p_Vid->pNextSlice->pic_parameter_set_id);
      MakePPSavailable (p_Vid, p_Vid->pNextPPS->pic_parameter_set_id, p_Vid->pNextPPS);
      p_Vid->pNextPPS->Valid=0;
    }

    //get the first slice from currentslice;
    assert(ppSliceList[p_Vid->iSliceNumOfCurrPic]);
    currSlice = ppSliceList[p_Vid->iSliceNumOfCurrPic];
    ppSliceList[p_Vid->iSliceNumOfCurrPic] = p_Vid->pNextSlice;
    p_Vid->pNextSlice = currSlice;
    assert(ppSliceList[p_Vid->iSliceNumOfCurrPic]->current_slice_nr == 0);
    
    currSlice = ppSliceList[p_Vid->iSliceNumOfCurrPic];

    UseParameterSet (currSlice);

    init_picture(p_Vid, currSlice, p_Inp);	//当读取到新一帧的第一个片时调用
    
    p_Vid->iSliceNumOfCurrPic++;
    current_header = SOS;
  }

	//从码流中将一个图片的所有条带信息读出,读到ppSliceList列表中
  while(current_header != SOP && current_header != EOS)	//= 3,4才执行
  {
    //no pending slices;
    assert(p_Vid->iSliceNumOfCurrPic < p_Vid->iNumOfSlicesAllocated);
    if(!ppSliceList[p_Vid->iSliceNumOfCurrPic]) //当前图片的条带列表
    {
      ppSliceList[p_Vid->iSliceNumOfCurrPic] = malloc_slice(p_Inp, p_Vid);
    }
    currSlice = ppSliceList[p_Vid->iSliceNumOfCurrPic];

    //p_Vid->currentSlice = currSlice;
    currSlice->p_Vid = p_Vid;
    currSlice->p_Inp = p_Inp;
    //currSlice->p_Dpb = p_Vid->p_Dpb_layer[0]; //set default value;
    currSlice->next_header = -8888;
    currSlice->num_dec_mb = 0;
    currSlice->coeff_ctr = -1;
    currSlice->pos       =  0;
    currSlice->is_reset_coeff = FALSE;
    currSlice->is_reset_coeff_cr = FALSE;

	//读取一个新的条带,每个条带打包成一个NALU,并处理slice_header
    current_header = read_new_slice(currSlice);	//大部分是SOP=2 最后是EOS	
    
    //init;
    currSlice->current_header = current_header;

    // error tracking of primary and redundant slices.
#if 0    
    Error_tracking(p_Vid, currSlice);
#endif
    // If primary and redundant are received and primary is correct, discard the redundant
    // else, primary slice will be replaced with redundant slice.
    //redundant_pic_cnt 对于那些属于基本编码图像的条带和条带数据分割块应等于0(非0冗余编码图像)
    if(currSlice->frame_num == p_Vid->previous_frame_num && currSlice->redundant_pic_cnt !=0
      && p_Vid->Is_primary_correct !=0 && current_header != EOS)
    {
      continue;
    }

    if((current_header != SOP && current_header != EOS) || (p_Vid->iSliceNumOfCurrPic == 0 && current_header == SOP))
    {
       currSlice->current_slice_nr = (short) p_Vid->iSliceNumOfCurrPic;
       p_Vid->dec_picture->max_slice_id = (short) imax(currSlice->current_slice_nr, p_Vid->dec_picture->max_slice_id);
       if(p_Vid->iSliceNumOfCurrPic >0)
       {
         //CopyPOC(*ppSliceList, currSlice);	//*ppSliceList -> currSlice
         ppSliceList[p_Vid->iSliceNumOfCurrPic-1]->end_mb_nr_plus1 = currSlice->start_mb_nr;
       }
	   
       p_Vid->iSliceNumOfCurrPic++;	//关键

	   //增加一个条带的分配 在alloc_video_params分配过一次
       if(p_Vid->iSliceNumOfCurrPic >= p_Vid->iNumOfSlicesAllocated)	
       {
         Slice **tmpSliceList = (Slice **)realloc(p_Vid->ppSliceList, (p_Vid->iNumOfSlicesAllocated+MAX_NUM_DECSLICES)*sizeof(Slice*));
         if(!tmpSliceList)
         {
           tmpSliceList = calloc((p_Vid->iNumOfSlicesAllocated+MAX_NUM_DECSLICES), sizeof(Slice*));
           memcpy(tmpSliceList, p_Vid->ppSliceList, p_Vid->iSliceNumOfCurrPic*sizeof(Slice*));
           //free;
           free(p_Vid->ppSliceList);
           ppSliceList = p_Vid->ppSliceList = tmpSliceList;
         }
         else
         {
           //assert(tmpSliceList == p_Vid->ppSliceList);
           ppSliceList = p_Vid->ppSliceList = tmpSliceList;
           memset(p_Vid->ppSliceList+p_Vid->iSliceNumOfCurrPic, 0, sizeof(Slice*)*MAX_NUM_DECSLICES);
         }
         p_Vid->iNumOfSlicesAllocated += MAX_NUM_DECSLICES;
       }

		p_Dec->nalu_pos_array_idx--;
       current_header = SOS;     //此处重新进入读取一个条带的过程(读取一个图片中得另外一个条带的过程)  
    }
    else
    {
      if(ppSliceList[p_Vid->iSliceNumOfCurrPic-1]->mb_aff_frame_flag)
       ppSliceList[p_Vid->iSliceNumOfCurrPic-1]->end_mb_nr_plus1 = p_Vid->FrameSizeInMbs/2;
      else
       ppSliceList[p_Vid->iSliceNumOfCurrPic-1]->end_mb_nr_plus1 = p_Vid->FrameSizeInMbs/(1+ppSliceList[p_Vid->iSliceNumOfCurrPic-1]->field_pic_flag);
	  
       p_Vid->newframe = 1;
       currSlice->current_slice_nr = 0;
	   
       //keep it in currentslice;
       ppSliceList[p_Vid->iSliceNumOfCurrPic] = p_Vid->pNextSlice;	//ppSliceList[p_Vid->iSliceNumOfCurrPic]地址变成0
       p_Vid->pNextSlice = currSlice; 	//p_Vid->pNextSlice指向下一个条带

	   p_Dec->nalu_pos_array_idx++;
    }
	
	#if 0
    copy_slice_info(currSlice, p_Vid->old_slice); //currSlice -> old_slice
	#endif
  }
	
  iRet = current_header;
  init_picture_decoding(p_Vid);		//需要修改

	//循环解码一帧中的所有条带
  for(iSliceNo=0; iSliceNo<p_Vid->iSliceNumOfCurrPic; iSliceNo++)
  {
    currSlice = ppSliceList[iSliceNo];
    current_header = currSlice->current_header;
    //p_Vid->currentSlice = currSlice;

    assert(current_header != EOS);
    assert(currSlice->current_slice_nr == iSliceNo);

    init_slice(p_Vid, currSlice);	//需要修改
	//printf("start3decode: idx:%2d, cur_rbsp_pos: %5d\n",p_Dec->nalu_pos_array_idx,p_Dec->nalu_pos_array[p_Dec->nalu_pos_array_idx]);
    decode_slice(currSlice, current_header);

    p_Vid->iNumOfSlicesDecoded++;
    p_Vid->num_dec_mb += currSlice->num_dec_mb;
    //p_Vid->erc_mvperMB += currSlice->erc_mvperMB;
  }
		
#if MVC_EXTENSION_ENABLE
  p_Vid->last_dec_view_id = p_Vid->dec_picture->view_id;
#endif
  //if(p_Vid->dec_picture->structure == FRAME)
    //p_Vid->last_dec_poc = p_Vid->dec_picture->frame_poc;
  //else if(p_Vid->dec_picture->structure == TOP_FIELD)
    //p_Vid->last_dec_poc = p_Vid->dec_picture->top_poc;
  //else if(p_Vid->dec_picture->structure == BOTTOM_FIELD)
    //p_Vid->last_dec_poc = p_Vid->dec_picture->bottom_poc;

#if 1	  
  exit_picture(p_Vid, &p_Vid->dec_picture);
  p_Vid->previous_frame_num = ppSliceList[0]->frame_num;
#endif 	

  return (iRet);	//返回current_header 最后一个宏块处理完后此处返回有问题,返回了2 SOP
}

/*!
 ************************************************************************
 * \brief
 *    Reads new slice from bit_stream_dec
 ************************************************************************
 */
int read_new_slice(Slice *currSlice)
{
  VideoParameters *p_Vid = currSlice->p_Vid;
  InputParameters *p_Inp = currSlice->p_Inp;

  NALU_t *nalu = p_Vid->nalu; 
  int current_header = 0;
  int BitsUsedByHeader;
  Bitstream *currStream = NULL;

  int slice_id_a, slice_id_b, slice_id_c;

  for (;;)	//处理完SPS,PPS,SEI后会break回到此处
  {
#if (MVC_EXTENSION_ENABLE)
    currSlice->svc_extension_flag = -1;
#endif
    //读取码流(test.264)获取NALU,并转化为RBSP,并返回了一个RBSP的长度(后续需处理的数据在nalu->buf中,长度是nalu->len)
    if (0 == read_next_nalu(p_Vid, nalu))  
      return EOS;

#if (MVC_EXTENSION_ENABLE)
    if(p_Inp->DecodeAllLayers == 1 && (nalu->nal_unit_type == NALU_TYPE_PREFIX || nalu->nal_unit_type == NALU_TYPE_SLC_EXT))
    {
      currStream = currSlice->partArr[0].bitstream;
      currStream->ei_flag = 0;
      currStream->frame_bitoffset = currStream->read_len = 0;
      fast_memcpy (currStream->streamBuffer, &nalu->buf[1], nalu->len-1);
      currStream->code_len = currStream->bitstream_length = RBSPtoSODB(currStream->streamBuffer, nalu->len-1);

      currSlice->svc_extension_flag = read_u_1 ("svc_extension_flag"        , currStream, p_Dec);

      if(currSlice->svc_extension_flag)
      {
        nal_unit_header_svc_extension();
      }
      else
      {
        nal_unit_header_mvc_extension(&currSlice->NaluHeaderMVCExt, currStream);
        currSlice->NaluHeaderMVCExt.iPrefixNALU = (nalu->nal_unit_type == NALU_TYPE_PREFIX);
      }

      if(nalu->nal_unit_type == NALU_TYPE_SLC_EXT)
      {        
        if(currSlice->svc_extension_flag)
        {
          //to be implemented for Annex G;
        }
        else 
        {
          nalu->nal_unit_type = NALU_TYPE_SLICE; //currSlice->NaluHeaderMVCExt.non_idr_flag==0? NALU_TYPE_IDR: NALU_TYPE_SLICE; 
        }
      }
    }
#endif

process_nalu:
    switch (nalu->nal_unit_type)
    {
    case NALU_TYPE_SLICE:
    case NALU_TYPE_IDR: //进入普通解码过程
		//slice都不使用(A/B/C)分区的过程
      if (p_Vid->recovery_point || nalu->nal_unit_type == NALU_TYPE_IDR)
      {
        if (p_Vid->recovery_point_found == 0)
        {
          if (nalu->nal_unit_type != NALU_TYPE_IDR)
          {
            printf("Warning: Decoding does not start with an IDR picture.\n");
            p_Vid->non_conforming_stream = 1;
          }
          else
            p_Vid->non_conforming_stream = 0;
        }
        p_Vid->recovery_point_found = 1;
      }

      if (p_Vid->recovery_point_found == 0)
        break;

      currSlice->idr_flag = (nalu->nal_unit_type == NALU_TYPE_IDR);
      currSlice->nal_reference_idc = nalu->nal_reference_idc;
      currSlice->dp_mode = PAR_DP_1;
      currSlice->max_part_nr = 1;
#if (MVC_EXTENSION_ENABLE)
      if (currSlice->svc_extension_flag != 0)  // = -1
      {
        currStream = currSlice->partArr[0].bitstream;
        currStream->ei_flag = 0;
        currStream->frame_bitoffset = currStream->read_len = 0;
        fast_memcpy (currStream->streamBuffer, &nalu->buf[1], nalu->len-1);	//将nalu除了头的数据拷贝给streamBuffer
        currStream->code_len = currStream->bitstream_length = RBSPtoSODB(currStream->streamBuffer, nalu->len-1);
      }
#else   
      currStream = currSlice->partArr[0].bitstream;
      currStream->ei_flag = 0;
      currStream->frame_bitoffset = currStream->read_len = 0;
      memcpy (currStream->streamBuffer, &nalu->buf[1], nalu->len-1);  //将nalu除了头的数据拷贝给streamBuffer
      currStream->code_len = currStream->bitstream_length = RBSPtoSODB(currStream->streamBuffer, nalu->len-1);
#endif

#if (MVC_EXTENSION_ENABLE)
      if(currSlice->svc_extension_flag == 0)
      {  //MVC:multiple view coding 用于支持3D视频和FTV多视点编码
        //if(is_MVC_profile(p_Vid->active_sps->profile_idc))
        //{
          currSlice->view_id = currSlice->NaluHeaderMVCExt.view_id;
          currSlice->inter_view_flag = currSlice->NaluHeaderMVCExt.inter_view_flag;
          currSlice->anchor_pic_flag = currSlice->NaluHeaderMVCExt.anchor_pic_flag;
        //}
      }
      else if(currSlice->svc_extension_flag == -1) //SVC(Scalable Video Coding可分级视频编码) and the normal AVC;
      {
        if(p_Vid->active_subset_sps == NULL)
        {
          currSlice->view_id = GetBaseViewId(p_Vid, &p_Vid->active_subset_sps);
          if(currSlice->NaluHeaderMVCExt.iPrefixNALU >0)
          {
            assert(currSlice->view_id == currSlice->NaluHeaderMVCExt.view_id);
            currSlice->inter_view_flag = currSlice->NaluHeaderMVCExt.inter_view_flag;
            currSlice->anchor_pic_flag = currSlice->NaluHeaderMVCExt.anchor_pic_flag;
          }
          else
          {
            currSlice->inter_view_flag = 1;
            currSlice->anchor_pic_flag = currSlice->idr_flag;
          }
        }
        else
        {
          assert(p_Vid->active_subset_sps->num_views_minus1 >=0);
          // prefix NALU available
          if(currSlice->NaluHeaderMVCExt.iPrefixNALU >0)
          {
            currSlice->view_id = currSlice->NaluHeaderMVCExt.view_id;
            currSlice->inter_view_flag = currSlice->NaluHeaderMVCExt.inter_view_flag;
            currSlice->anchor_pic_flag = currSlice->NaluHeaderMVCExt.anchor_pic_flag;
          }
          else
          { //no prefix NALU;
            currSlice->view_id = p_Vid->active_subset_sps->view_id[0];
            currSlice->inter_view_flag = 1;
            currSlice->anchor_pic_flag = currSlice->idr_flag;
          }
        }
      }
     currSlice->layer_id = currSlice->view_id = GetVOIdx( p_Vid, currSlice->view_id );
      /*
      if(currSlice->view_id >=0)
      {
        currSlice->p_Dpb = p_Vid->p_Dpb_layer[currSlice->view_id];
        currSlice->imgtype = currSlice->p_Dpb->imgtype;
      }
      */
#endif

      // Some syntax of the Slice Header depends on the parameter set, which depends on
      // the parameter set ID of the SLice header.  Hence, read the pic_parameter_set_id
      // of the slice header first, then setup the active parameter sets, and then read
      // the rest of the slice header
      //条带头语法
      BitsUsedByHeader = FirstPartOfSliceHeader(currSlice);  //读取条带头的三个语法元素(first_mb_in_slice~pps_id)
      UseParameterSet (currSlice);	//设置参数
      currSlice->active_sps = p_Vid->active_sps;
      currSlice->active_pps = p_Vid->active_pps;
      currSlice->Transform8x8Mode = p_Vid->active_pps->transform_8x8_mode_flag;
      currSlice->chroma444_not_separate = (p_Vid->active_sps->chroma_format_idc==YUV444)&&((p_Vid->separate_colour_plane_flag == 0));

      BitsUsedByHeader = RestOfSliceHeader (currSlice);	//BitsUsedByHeader包括了slice_header语法元素的位数

			/* record slice_header used bits *///从p_Dec->cur_nal_start_pos+1+BitsUsedByHeader/8开始为当前条带的第一个宏块的数据
#if (MVC_EXTENSION_ENABLE)
      if(currSlice->view_id >=0)
      {
        //currSlice->p_Dpb = p_Vid->p_Dpb_layer[currSlice->view_id];	//解码图像缓存
      }
#endif

      assign_quant_params (currSlice);	//分配量化参数     

      // if primary slice is replaced with redundant slice, set the correct image type
      if(currSlice->redundant_pic_cnt && p_Vid->Is_primary_correct==0 && p_Vid->Is_redundant_correct)
      {
        p_Vid->dec_picture->slice_type = p_Vid->type;
      }

      if(is_new_picture(p_Vid->dec_picture, currSlice, p_Vid->old_slice))
      {
        if(p_Vid->iSliceNumOfCurrPic==0)	//如果是一个新的图像且是第一个条带,则初始化一个图片
          init_picture(p_Vid, currSlice, p_Inp);

        current_header = SOP;
        //check zero_byte if it is also the first NAL unit in the access unit
        CheckZeroByteVCL(p_Vid, nalu);
      }
      else
        current_header = SOS;

		/****根据slice type设置预测模式(read mv info from NAL)*/
		//包含slice_data()语法
      setup_slice_methods(currSlice);	

      // From here on, p_Vid->active_sps, p_Vid->active_pps and the slice header are valid
      if (currSlice->mb_aff_frame_flag)
        currSlice->current_mb_nr = currSlice->start_mb_nr << 1;
      else
        currSlice->current_mb_nr = currSlice->start_mb_nr;

      if (p_Vid->active_pps->entropy_coding_mode_flag)
      {
        int ByteStartPosition = currStream->frame_bitoffset/8;
        if (currStream->frame_bitoffset%8 != 0)
        {
          ++ByteStartPosition;
        }
        arideco_start_decoding (&currSlice->partArr[0].de_cabac, currStream->streamBuffer, ByteStartPosition, &currStream->read_len);
      }

      //FreeNALU(nalu);
      p_Vid->recovery_point = 0;

      return current_header;
      break;
    case NALU_TYPE_DPA:
	  //处理顺序:A->B->C
      if (p_Vid->recovery_point_found == 0)
        break;

      // read DP_A
      currSlice->dpB_NotPresent =1; 
      currSlice->dpC_NotPresent =1; 

      currSlice->idr_flag          = FALSE;
      currSlice->nal_reference_idc = nalu->nal_reference_idc;
      currSlice->dp_mode     = PAR_DP_3;
      currSlice->max_part_nr = 3;
      currSlice->ei_flag     = 0;
#if MVC_EXTENSION_ENABLE
      //currSlice->p_Dpb = p_Vid->p_Dpb_layer[0];
#endif
      currStream             = currSlice->partArr[0].bitstream;
      currStream->ei_flag    = 0;
      currStream->frame_bitoffset = currStream->read_len = 0;
      memcpy (currStream->streamBuffer, &nalu->buf[1], nalu->len-1);	//将nalu除了头的数据拷贝给streamBuffer
      currStream->code_len = currStream->bitstream_length = RBSPtoSODB(currStream->streamBuffer, nalu->len-1);
#if MVC_EXTENSION_ENABLE
      currSlice->view_id = GetBaseViewId(p_Vid, &p_Vid->active_subset_sps);
      currSlice->inter_view_flag = 1;
      currSlice->layer_id = currSlice->view_id = GetVOIdx( p_Vid, currSlice->view_id );
      currSlice->anchor_pic_flag = currSlice->idr_flag;
#endif

      BitsUsedByHeader = FirstPartOfSliceHeader(currSlice);
      UseParameterSet (currSlice);
      currSlice->active_sps = p_Vid->active_sps;
      currSlice->active_pps = p_Vid->active_pps;
      currSlice->Transform8x8Mode = p_Vid->active_pps->transform_8x8_mode_flag;
      currSlice->chroma444_not_separate = (p_Vid->active_sps->chroma_format_idc==YUV444)&&((p_Vid->separate_colour_plane_flag == 0));

      BitsUsedByHeader += RestOfSliceHeader (currSlice);
#if MVC_EXTENSION_ENABLE
      //currSlice->p_Dpb = p_Vid->p_Dpb_layer[currSlice->view_id];
#endif

      assign_quant_params (currSlice);        

      if(is_new_picture(p_Vid->dec_picture, currSlice, p_Vid->old_slice))
      {
        if(p_Vid->iSliceNumOfCurrPic==0)
          init_picture(p_Vid, currSlice, p_Inp);

        current_header = SOP;
        //check zero_byte if it is also the first NAL unit in the access unit
        CheckZeroByteVCL(p_Vid, nalu);
      }
      else
        current_header = SOS;

      setup_slice_methods(currSlice);

      // From here on, p_Vid->active_sps, p_Vid->active_pps and the slice header are valid
      if (currSlice->mb_aff_frame_flag)
        currSlice->current_mb_nr = currSlice->start_mb_nr << 1;
      else
        currSlice->current_mb_nr = currSlice->start_mb_nr;

      // Now I need to read the slice ID, which depends on the value of
      // redundant_pic_cnt_present_flag

      slice_id_a  = read_ue_v("NALU: DP_A slice_id", currStream, p_Dec);	//读取slice_id

      if (p_Vid->active_pps->entropy_coding_mode_flag)
        error ("received data partition with CABAC, this is not allowed", 500);

      // continue with reading next DP(Data Partition)
      if (0 == read_next_nalu(p_Vid, nalu))
        return current_header;

      if ( NALU_TYPE_DPB == nalu->nal_unit_type)
      {
        // we got a DPB
        currStream             = currSlice->partArr[1].bitstream;
        currStream->ei_flag    = 0;
        currStream->frame_bitoffset = currStream->read_len = 0;

        memcpy (currStream->streamBuffer, &nalu->buf[1], nalu->len-1);
        currStream->code_len = currStream->bitstream_length = RBSPtoSODB(currStream->streamBuffer, nalu->len-1);

        slice_id_b  = read_ue_v("NALU: DP_B slice_id", currStream, p_Dec);

        currSlice->dpB_NotPresent = 0; 

        if ((slice_id_b != slice_id_a) || (nalu->lost_packets))
        {
          printf ("Waning: got a data partition B which does not match DP_A (DP loss!)\n");
          currSlice->dpB_NotPresent =1; 
          currSlice->dpC_NotPresent =1; 
        }
        else
        {
          if (p_Vid->active_pps->redundant_pic_cnt_present_flag)
            read_ue_v("NALU: DP_B redundant_pic_cnt", currStream, p_Dec);

          // we're finished with DP_B, so let's continue with next DP
          if (0 == read_next_nalu(p_Vid, nalu))
            return current_header;
        }
      }
      else
      {
        currSlice->dpB_NotPresent =1; 
      }

      // check if we got DP_C
      if ( NALU_TYPE_DPC == nalu->nal_unit_type)
      {
        currStream             = currSlice->partArr[2].bitstream;
        currStream->ei_flag    = 0;
        currStream->frame_bitoffset = currStream->read_len = 0;

        memcpy (currStream->streamBuffer, &nalu->buf[1], nalu->len-1);
        currStream->code_len = currStream->bitstream_length = RBSPtoSODB(currStream->streamBuffer, nalu->len-1);

        currSlice->dpC_NotPresent = 0;

        slice_id_c  = read_ue_v("NALU: DP_C slice_id", currStream, p_Dec);
        if ((slice_id_c != slice_id_a)|| (nalu->lost_packets))
        {
          printf ("Warning: got a data partition C which does not match DP_A(DP loss!)\n");
          //currSlice->dpB_NotPresent =1;
          currSlice->dpC_NotPresent =1;
        }

        if (p_Vid->active_pps->redundant_pic_cnt_present_flag)
          read_ue_v("NALU:SLICE_C redudand_pic_cnt", currStream, p_Dec);
      }
      else
      {
        currSlice->dpC_NotPresent =1;
      }

      // check if we read anything else than the expected partitions
      if ((nalu->nal_unit_type != NALU_TYPE_DPB) && (nalu->nal_unit_type != NALU_TYPE_DPC))
      {
        // we have a NALI that we can't process here, so restart processing
        goto process_nalu;
        // yes, "goto" should not be used, but it's really the best way here before we restructure the decoding loop
        // (which should be taken care of anyway)
      }

      //FreeNALU(nalu);
      return current_header;
      break;
    case NALU_TYPE_DPB:
      if (p_Inp->silent == FALSE)
      {
        printf ("found data partition B without matching DP A, discarding\n");
      }
      break;
    case NALU_TYPE_DPC:
      if (p_Inp->silent == FALSE)
      {
        printf ("found data partition C without matching DP A, discarding\n");
      }
      break;
    case NALU_TYPE_SEI:
		p_Dec->nalu_pos_array_idx++;
      InterpretSEIMessage(nalu->buf,nalu->len,p_Vid, currSlice);
      break;
    case NALU_TYPE_PPS:
		p_Dec->nalu_pos_array_idx++;
      ProcessPPS(p_Vid, nalu);
      break;
    case NALU_TYPE_SPS:
		p_Dec->nalu_pos_array_idx++;
      ProcessSPS(p_Vid, nalu);
      break;
    case NALU_TYPE_AUD:
      break;
    case NALU_TYPE_EOSEQ:
      break;
    case NALU_TYPE_EOSTREAM:
      break;
    case NALU_TYPE_FILL:
      if (p_Inp->silent == FALSE)
      {
        //printf ("read_new_slice: Found NALU_TYPE_FILL, len %d\n", (int) nalu->len);
        //printf ("Skipping these filling bits, proceeding w/ next NALU\n");
      }
      break;
#if (MVC_EXTENSION_ENABLE)
    case NALU_TYPE_VDRD:
			break;
    case NALU_TYPE_PREFIX:
      if(currSlice->svc_extension_flag==1)
        prefix_nal_unit_svc();	//empty function
      break;
    case NALU_TYPE_SUB_SPS:
      if (p_Inp->DecodeAllLayers== 1)
      {
        ProcessSubsetSPS(p_Vid, nalu);
      }
      else
      {
        if (p_Inp->silent == FALSE)
          printf ("Found Subsequence SPS NALU. Ignoring.\n");
      }
      break;
    case NALU_TYPE_SLC_EXT:
      if (p_Inp->DecodeAllLayers == 0 &&  (p_Inp->silent == FALSE))
        printf ("Found SVC extension NALU (%d). Ignoring.\n", (int) nalu->nal_unit_type);
      break;
#endif
    default:
      {
        if (p_Inp->silent == FALSE)
          printf ("Found NALU type %d, len %d undefined, ignore NALU, moving on\n", (int) nalu->nal_unit_type, (int) nalu->len);
      }
      break;
    }
  }
}


/*!
 ************************************************************************
 * \brief
 *    finish decoding of a picture, conceal errors and store it
 *    into the DPB
 ************************************************************************
 */
void exit_picture(VideoParameters *p_Vid, StorablePicture **dec_picture)
{
  InputParameters *p_Inp = p_Vid->p_Inp;
  //SNRParameters   *snr   = p_Vid->snr;
  char yuv_types[4][6]= {"4:0:0","4:2:0","4:2:2","4:4:4"};

  int structure, slice_type, refpic, qp, pic_num, chroma_format_idc, is_idr;

  int64 tmp_time;                   // time used by decoding the last frame
  char   yuvFormat[10];

  // return if the last picture has already been finished
  if (*dec_picture==NULL || (p_Vid->num_dec_mb != p_Vid->PicSizeInMbs && (p_Vid->yuv_format != YUV444 || !p_Vid->separate_colour_plane_flag)))
  {
    return;
  }
#if 0
  //if(p_Vid->bDeblockEnable & (1<<(*dec_picture)->used_for_reference))
  if(0x03 & (1<<(*dec_picture)->used_for_reference))
  {
    //deblocking for frame or field
    if( (p_Vid->separate_colour_plane_flag != 0) )
    {
      int nplane;
      int colour_plane_id = p_Vid->ppSliceList[0]->colour_plane_id;
      for( nplane=0; nplane<MAX_PLANE; ++nplane )
      {
        p_Vid->ppSliceList[0]->colour_plane_id = nplane;
        change_plane_JV( p_Vid, nplane, NULL );
        //DeblockPicture( p_Vid, *dec_picture );
      }
      p_Vid->ppSliceList[0]->colour_plane_id = colour_plane_id;
      make_frame_picture_JV(p_Vid);
    }
    else
    {
      //DeblockPicture( p_Vid, *dec_picture );
    }
  }
  else
  {
    if( (p_Vid->separate_colour_plane_flag != 0) )
    {
      make_frame_picture_JV(p_Vid);
    }
  }
#endif
  //if ((*dec_picture)->mb_aff_frame_flag)
    //MbAffPostProc(p_Vid);

  if (p_Vid->structure == FRAME)         // buffer mgt. for frame mode
    frame_postprocessing(p_Vid);
  else
    field_postprocessing(p_Vid);   // reset all interlaced variables

#if 0
  structure  = (*dec_picture)->structure;
  slice_type = (*dec_picture)->slice_type;
  //frame_poc  = (*dec_picture)->frame_poc;  
  refpic     = (*dec_picture)->used_for_reference;
  qp         = (*dec_picture)->qp;
  pic_num    = (*dec_picture)->pic_num;
  is_idr     = (*dec_picture)->idr_flag;

  chroma_format_idc = (*dec_picture)->chroma_format_idc;

  *dec_picture=NULL;
#endif 

  if (p_Vid->last_has_mmco_5)
  {
    p_Vid->pre_frame_num = 0;
  }

  if (p_Inp->silent == FALSE)
  {
    if (structure==TOP_FIELD || structure==FRAME)
    {
      if(slice_type == I_SLICE && is_idr) // IDR picture
        strcpy(p_Vid->cslice_type,"IDR");
      else if(slice_type == I_SLICE) // I picture
        strcpy(p_Vid->cslice_type," I ");
      else if(slice_type == P_SLICE) // P pictures
        strcpy(p_Vid->cslice_type," P ");
      else if(slice_type == SP_SLICE) // SP pictures
        strcpy(p_Vid->cslice_type,"SP ");
      else if (slice_type == SI_SLICE)
        strcpy(p_Vid->cslice_type,"SI ");
      else if(refpic) // stored B pictures
        strcpy(p_Vid->cslice_type," B ");
      else // B pictures
        strcpy(p_Vid->cslice_type," b ");

      if (structure==FRAME)
      {
        strncat(p_Vid->cslice_type,")       ",8-strlen(p_Vid->cslice_type));
      }
    }
    else if (structure==BOTTOM_FIELD)
    {
      if(slice_type == I_SLICE && is_idr) // IDR picture
        strncat(p_Vid->cslice_type,"|IDR)",8-strlen(p_Vid->cslice_type));
      else if(slice_type == I_SLICE) // I picture
        strncat(p_Vid->cslice_type,"| I )",8-strlen(p_Vid->cslice_type));
      else if(slice_type == P_SLICE) // P pictures
        strncat(p_Vid->cslice_type,"| P )",8-strlen(p_Vid->cslice_type));
      else if(slice_type == SP_SLICE) // SP pictures
        strncat(p_Vid->cslice_type,"|SP )",8-strlen(p_Vid->cslice_type));
      else if (slice_type == SI_SLICE)
        strncat(p_Vid->cslice_type,"|SI )",8-strlen(p_Vid->cslice_type));
      else if(refpic) // stored B pictures
        strncat(p_Vid->cslice_type,"| B )",8-strlen(p_Vid->cslice_type));
      else // B pictures
        strncat(p_Vid->cslice_type,"| b )",8-strlen(p_Vid->cslice_type));   
    }
  }

  if ((structure==FRAME)||structure==BOTTOM_FIELD)
  {
    gettime (&(p_Vid->end_time));              // end time

    tmp_time  = timediff(&(p_Vid->start_time), &(p_Vid->end_time));
    p_Vid->tot_time += tmp_time;
    tmp_time  = timenorm(tmp_time);
    sprintf(yuvFormat,"%s", yuv_types[chroma_format_idc]);

#if 1
    if (p_Inp->silent == FALSE)
    {
      //SNRParameters   *snr = p_Vid->snr;
      //if (p_Vid->p_ref != -1)
        //fprintf(stdout,"%05d(%s%5d %5d %5d %8.4f %8.4f %8.4f  %s %7d\n",
        //p_Vid->frame_no, p_Vid->cslice_type, frame_poc, pic_num, qp, snr->snr[0], snr->snr[1], snr->snr[2], yuvFormat, (int) tmp_time);
      //else
      //此处输出每帧的信息
        fprintf(stdout,"%05d(%s%5d %5d                             %s %7d\n",
        p_Vid->frame_no, p_Vid->cslice_type, pic_num, qp, yuvFormat, (int)tmp_time);
    }
    //else
      //fprintf(stdout,"Completed Decoding frame %05d.\r",snr->frame_ctr);
#endif

    fflush(stdout);
	

    if(slice_type == I_SLICE || slice_type == SI_SLICE || slice_type == P_SLICE || refpic)   // I or P pictures
    {
#if (MVC_EXTENSION_ENABLE)
      if((p_Vid->ppSliceList[0])->view_id!=0)
#endif
        ++(p_Vid->number);
    }
    else
      ++(p_Vid->Bframe_ctr);    // B pictures
    //++(snr->frame_ctr);

#if (MVC_EXTENSION_ENABLE)
    if ((p_Vid->ppSliceList[0])->view_id != 0)
#endif
      ++(p_Vid->g_nFrame);   
  }

  //p_Vid->currentSlice->current_mb_nr = -4712;   // impossible value for debugging, StW
  //p_Vid->currentSlice->current_slice_nr = 0;
}

/*!
 ************************************************************************
 * \brief
 *    set defaults for old_slice
 *    NAL unit of a picture"
 ************************************************************************
 */
void init_old_slice(OldSliceParams *p_old_slice)
{
  p_old_slice->field_pic_flag = 0;
  p_old_slice->pps_id         = INT_MAX;
  p_old_slice->frame_num      = INT_MAX;
  p_old_slice->nal_ref_idc    = INT_MAX;
  p_old_slice->idr_flag       = FALSE;

  p_old_slice->pic_oder_cnt_lsb          = UINT_MAX;
  p_old_slice->delta_pic_oder_cnt_bottom = INT_MAX;

  p_old_slice->delta_pic_order_cnt[0] = INT_MAX;
  p_old_slice->delta_pic_order_cnt[1] = INT_MAX;
}

//copy currSlice to p_old_slice
void copy_slice_info(Slice *currSlice, OldSliceParams *p_old_slice)
{
  VideoParameters *p_Vid = currSlice->p_Vid;

  p_old_slice->pps_id         = currSlice->pic_parameter_set_id;
  p_old_slice->frame_num      = currSlice->frame_num; //p_Vid->frame_num;
  p_old_slice->field_pic_flag = currSlice->field_pic_flag; //p_Vid->field_pic_flag;

  if(currSlice->field_pic_flag)
  {
    p_old_slice->bottom_field_flag = currSlice->bottom_field_flag;
  }

  p_old_slice->nal_ref_idc = currSlice->nal_reference_idc;
  p_old_slice->idr_flag    = (byte) currSlice->idr_flag;

  if (currSlice->idr_flag)
  {
    p_old_slice->idr_pic_id = currSlice->idr_pic_id;
  }

  if (p_Vid->active_sps->pic_order_cnt_type == 0)
  {
    p_old_slice->pic_oder_cnt_lsb          = currSlice->pic_order_cnt_lsb;
    p_old_slice->delta_pic_oder_cnt_bottom = currSlice->delta_pic_order_cnt_bottom;
  }

  if (p_Vid->active_sps->pic_order_cnt_type == 1)
  {
    p_old_slice->delta_pic_order_cnt[0] = currSlice->delta_pic_order_cnt[0];
    p_old_slice->delta_pic_order_cnt[1] = currSlice->delta_pic_order_cnt[1];
  }
#if (MVC_EXTENSION_ENABLE)
  p_old_slice->view_id = currSlice->view_id;
  p_old_slice->inter_view_flag = currSlice->inter_view_flag; 
  p_old_slice->anchor_pic_flag = currSlice->anchor_pic_flag;
#endif
  p_old_slice->layer_id = currSlice->layer_id;
}

/*!
 ************************************************************************
 * \brief
 *    detect if current slice is "first VCL NAL unit of a picture"
 ************************************************************************
 */
int is_new_picture(StorablePicture *dec_picture, Slice *currSlice, OldSliceParams *p_old_slice)
{
  VideoParameters *p_Vid = currSlice->p_Vid;

  int result=0;

  result |= (NULL==dec_picture);

  result |= (p_old_slice->pps_id != currSlice->pic_parameter_set_id);

  result |= (p_old_slice->frame_num != currSlice->frame_num);

  result |= (p_old_slice->field_pic_flag != currSlice->field_pic_flag);

  if(currSlice->field_pic_flag && p_old_slice->field_pic_flag)
  {
    result |= (p_old_slice->bottom_field_flag != currSlice->bottom_field_flag);
  }

  result |= (p_old_slice->nal_ref_idc != currSlice->nal_reference_idc) && ((p_old_slice->nal_ref_idc == 0) || (currSlice->nal_reference_idc == 0));
  result |= (p_old_slice->idr_flag    != currSlice->idr_flag);

  if (currSlice->idr_flag && p_old_slice->idr_flag)
  {
    result |= (p_old_slice->idr_pic_id != currSlice->idr_pic_id);
  }

  if (p_Vid->active_sps->pic_order_cnt_type == 0)
  {
    result |= (p_old_slice->pic_oder_cnt_lsb          != currSlice->pic_order_cnt_lsb);
    if( p_Vid->active_pps->bottom_field_pic_order_in_frame_present_flag  ==  1 &&  !currSlice->field_pic_flag )
    {
      result |= (p_old_slice->delta_pic_oder_cnt_bottom != currSlice->delta_pic_order_cnt_bottom);
    }
  }

  if (p_Vid->active_sps->pic_order_cnt_type == 1)
  {
    if (!p_Vid->active_sps->delta_pic_order_always_zero_flag)
    {
      result |= (p_old_slice->delta_pic_order_cnt[0] != currSlice->delta_pic_order_cnt[0]);
      if( p_Vid->active_pps->bottom_field_pic_order_in_frame_present_flag  ==  1 &&  !currSlice->field_pic_flag )
      {
        result |= (p_old_slice->delta_pic_order_cnt[1] != currSlice->delta_pic_order_cnt[1]);
      }
    }
  }

#if (MVC_EXTENSION_ENABLE)
  result |= (currSlice->view_id != p_old_slice->view_id);
  result |= (currSlice->inter_view_flag != p_old_slice->inter_view_flag);
  result |= (currSlice->anchor_pic_flag != p_old_slice->anchor_pic_flag);
#endif
  result |= (currSlice->layer_id != p_old_slice->layer_id);
  return result;
}

/*!
 ************************************************************************
 * \brief
 *    Prepare field and frame buffer after frame decoding
 ************************************************************************
 */
void frame_postprocessing(VideoParameters *p_Vid)
{
}

/*!
 ************************************************************************
 * \brief
 *    Prepare field and frame buffer after field decoding
 ************************************************************************
 */
void field_postprocessing(VideoParameters *p_Vid)
{
  p_Vid->number /= 2;
}

/*!
 ************************************************************************
 * \brief
 *    decodes one slice
 ************************************************************************
 */
void decode_one_slice(Slice *currSlice)
{
  VideoParameters *p_Vid = currSlice->p_Vid;
  Boolean end_of_slice = FALSE;
  Macroblock *currMB = NULL;
  currSlice->cod_counter=-1;

  if( (p_Vid->separate_colour_plane_flag != 0) )
  {
    //change_plane_JV( p_Vid, currSlice->colour_plane_id, currSlice );
  }
  else
  {
    currSlice->mb_data = p_Vid->mb_data;
    currSlice->dec_picture = p_Vid->dec_picture;
    currSlice->siblock = p_Vid->siblock;
    currSlice->ipredmode = p_Vid->ipredmode;
    currSlice->intra_block = p_Vid->intra_block;
  }

#if 0
  if (currSlice->slice_type == B_SLICE)
  {  	
    //为直接预测模式做一些准备工作:获取co_located图像,计算mv_scale
    //compute_colocated(currSlice, currSlice->listX); 	
  }
		
  if (currSlice->slice_type != I_SLICE && currSlice->slice_type != SI_SLICE)
    init_cur_imgy(currSlice,p_Vid); 
  //reset_ec_flags(p_Vid);
#endif

  //宏块解码循环
  while (end_of_slice == FALSE) // loop over macroblocks
  {

#if TRACE
    fprintf(p_Dec->p_trace,"\n*********** MB: %i Slice: %i Type %d **********\n", currSlice->current_mb_nr, currSlice->current_slice_nr, currSlice->slice_type);
#endif

    // Initializes the current macroblock
    //计算宏块，块，像素的坐标；宏块结构语法元素的初始化；
    //相邻块的可用性；以及滤波参数

    start_macroblock(currSlice, &currMB);
    // Get the syntax elements from the NAL
    //熵解码：包括解出宏块类型、预测模式、MVD、CBP、
    //残差（包括反量化操作）等
    //read_one_macroblock_i_slice_cabac read_one_macroblock_i_slice_cavlc 
    currSlice->read_one_macroblock(currMB);   

#if 0
   //反变换及运动补偿：反量化反变换、运动补偿、像素重构等
   decode_one_macroblock(currMB, currSlice->dec_picture);

    if(currSlice->mb_aff_frame_flag && currMB->mb_field)
    {
      currSlice->num_ref_idx_active[LIST_0] >>= 1;
      currSlice->num_ref_idx_active[LIST_1] >>= 1;
    }

#if (DISABLE_ERC == 0)
    //写入各个8*8块的预测模式及运动向量到错误隐藏变量中
    ercWriteMBMODEandMV(currMB);
#endif
#endif
	
	//每解完一个宏块，调用宏块解码后处理函数,判断是否解码完一个条带
    end_of_slice = exit_macroblock(currSlice, (!currSlice->mb_aff_frame_flag|| currSlice->current_mb_nr%2));
  }  //宏块循环
}

#if (MVC_EXTENSION_ENABLE)
int GetVOIdx(VideoParameters *p_Vid, int iViewId)
{
  int iVOIdx = -1;
  int *piViewIdMap;
  if(p_Vid->active_subset_sps)
  {
    piViewIdMap = p_Vid->active_subset_sps->view_id;
    for(iVOIdx = p_Vid->active_subset_sps->num_views_minus1; iVOIdx>=0; iVOIdx--)
      if(piViewIdMap[iVOIdx] == iViewId)
        break;
  }
  else
  {
    subset_seq_parameter_set_rbsp_t *curr_subset_sps;
    int i;

    curr_subset_sps = p_Vid->SubsetSeqParSet;
    for(i=0; i<MAXSPS; i++)
    {
      if(curr_subset_sps->num_views_minus1>=0 && curr_subset_sps->sps.Valid)
      {
        break;
      }
      curr_subset_sps++;
    }

    if( i < MAXSPS )
    {
      p_Vid->active_subset_sps = curr_subset_sps;

      piViewIdMap = p_Vid->active_subset_sps->view_id;
      for(iVOIdx = p_Vid->active_subset_sps->num_views_minus1; iVOIdx>=0; iVOIdx--)
        if(piViewIdMap[iVOIdx] == iViewId)
          break;

      return iVOIdx;
    }
    else
    {
      iVOIdx = 0;
    }
  }

  return iVOIdx;
}

int GetViewIdx(VideoParameters *p_Vid, int iVOIdx)
{
  int iViewIdx = -1;
  int *piViewIdMap;

  if( p_Vid->active_subset_sps )
  {
    assert( p_Vid->active_subset_sps->num_views_minus1 >= iVOIdx && iVOIdx >= 0 );
    piViewIdMap = p_Vid->active_subset_sps->view_id;
    iViewIdx = piViewIdMap[iVOIdx];    
  }

  return iViewIdx;
}
int get_maxViewIdx (VideoParameters *p_Vid, int view_id, int anchor_pic_flag, int listidx)
{
  int VOIdx;
  int maxViewIdx = 0;

  VOIdx = view_id; 
  if(VOIdx >= 0)
  {
    if(anchor_pic_flag)
      maxViewIdx = listidx? p_Vid->active_subset_sps->num_anchor_refs_l1[VOIdx] : p_Vid->active_subset_sps->num_anchor_refs_l0[VOIdx];
    else
      maxViewIdx = listidx? p_Vid->active_subset_sps->num_non_anchor_refs_l1[VOIdx] : p_Vid->active_subset_sps->num_non_anchor_refs_l0[VOIdx];
  }

  return maxViewIdx;
}
#endif
