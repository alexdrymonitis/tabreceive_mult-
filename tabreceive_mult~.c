/*********************************************************************************
* (c) 2017 Alexandros Drymonitis                                                 *
* (c) 2017 IOhannes m zmölnig                                                    *
* Based on code (c) 1997-1999 Miller Puckette                                    *
* All rights reserved.                                                           *
*                                                                                *
* Redistribution and use in source and binary forms, with or without             *
* modification, are permitted provided that the following conditions are met:    *
*                                                                                *
* 1. Redistributions of source code must retain the above copyright notice, this *
*    list of conditions and the following disclaimer.                            *
*                                                                                *
* 2. Redistributions in binary form must reproduce the above copyright notice,   *
*    this list of conditions and the following disclaimer in the documentation   *
*    and/or other materials provided with the distribution.                      *
*                                                                                *
* 3. Neither the name of the copyright holder nor the names of its               *
*    contributors may be used to endorse or promote products derived from        *
*    this software without specific prior written permission.                    *
*                                                                                *
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"    *
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE      *
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE *
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE   *
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL     *
* DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR     *
* SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER     *
* CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,  *
* OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE  *
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.           *
**********************************************************************************/

#include "m_pd.h"
/* include string.h to use memmove */
#include <string.h>

static t_class *tabreceive_mult_class;

/* define a data structure for each table added */
typedef struct _tabledata
{
  /* table data */
  t_word *x_vec;
  t_symbol *x_arrayname;
  int x_npoints;
  /* ramp data */
  t_sample x_target; /* target value of ramp */
  t_sample x_value; /* current value of ramp at block-borders */
  t_sample x_biginc;
  t_sample x_inc;
  t_float	x_time_ms;
  int x_ticksleft;
  int x_retarget;
} t_tabledata;

/* define the main data structure of the object */
typedef struct _tabreceive_mult
{
  t_object x_obj;
  t_tabledata *x_tables;
  t_float x_1overn;
  t_float x_dspticktomsec;
  t_int x_numtabs;
  int x_whichtab;
  /* set an int to test whether the object has just been created
  this way, if the object is in a patch with an array it listens to
  on load the array won't be found. this int protects from deleting the array
  when the object is created */
  int x_init;
} t_tabreceive_mult;

static t_int *tabreceive_mult_perform(t_int *w)
{
  t_tabreceive_mult *x = (t_tabreceive_mult *)(w[1]);
  t_sample *out = (t_sample *)(w[2]);
  int n = w[3];
  int i, j;

  /* first zero the output buffer */
  for (i = 0; i < n; i++) {
    out[i] = 0;
  }

  /* then run through each table and accumulate their values to the output */
  for (i = 0; i < x->x_numtabs; i++) {
    int vecsize;
    t_word *vec = x->x_tables[i].x_vec;

    /* check if the table doesn't exist
    use 0 == vec instead of vec == 0 because it's less error prone
    in case you forget to type two equal signs */
    if (0 == vec) continue;

    /* set vecsize to n and then check it against the size of each array */
    vecsize = n;
    if (vecsize > x->x_tables[i].x_npoints) {
      vecsize = x->x_tables[i].x_npoints;
    }

    /* set the values for the ramps */
    if (PD_BIGORSMALL(x->x_tables[i].x_value)) {
      x->x_tables[i].x_value = 0;
    }

    if (x->x_tables[i].x_retarget) {
      int nticks = (int)x->x_tables[i].x_time_ms * x->x_dspticktomsec;
      if (!nticks) nticks = 1;
      x->x_tables[i].x_ticksleft = nticks;
      x->x_tables[i].x_biginc = (x->x_tables[i].x_target - x->x_tables[i].x_value)/(t_float)nticks;
      x->x_tables[i].x_inc = x->x_1overn * x->x_tables[i].x_biginc;
      x->x_tables[i].x_retarget = 0;
    }
    if (x->x_tables[i].x_ticksleft) {
      t_sample f = x->x_tables[i].x_value;
      /* set the output values inside these tests */
      for (j = 0; j < vecsize; j++) {
        out[j] += vec[j].w_float * f, f += x->x_tables[i].x_inc;
      }
      x->x_tables[i].x_value += x->x_tables[i].x_biginc;
      x->x_tables[i].x_ticksleft--;
    }
    else {
      t_sample g = x->x_tables[i].x_value = x->x_tables[i].x_target;
      for (j = 0; j < vecsize; j++) {
        out[j] += vec[j].w_float * g;
      }
    }
  }
  return (w+4);
}

static void tabreceive_mult_setramp(t_tabreceive_mult *x, t_sample target, t_float time_ms)
{
  if (time_ms < 0.0f) time_ms = 0.0f;
  x->x_tables[x->x_whichtab].x_time_ms = time_ms;
  /* copied and modified from [line~]'s code */
  if (x->x_tables[x->x_whichtab].x_time_ms <= 0) {
    x->x_tables[x->x_whichtab].x_target = x->x_tables[x->x_whichtab].x_value = target;
    x->x_tables[x->x_whichtab].x_ticksleft = x->x_tables[x->x_whichtab].x_retarget = 0;
  }
  else {
    x->x_tables[x->x_whichtab].x_target = target;
    x->x_tables[x->x_whichtab].x_retarget = 1;
  }
}

static void tabreceive_mult_addtable(t_tabreceive_mult *x, t_symbol *s)
{
  int i;
  /* check if the array already exists in the list */
  for(i = 0; i < x->x_numtabs; i++) {
    if (x->x_tables[i].x_arrayname == s) {
      x->x_whichtab = i;
      return;
    }
  }

  /* if the array doesn't exist, reallocate memory for the new array */
  x->x_tables = (t_tabledata*)resizebytes(x->x_tables,
    sizeof(t_tabledata) * x->x_numtabs,
    sizeof(t_tabledata) * (x->x_numtabs+1));

  /* once the extra memory is reallocated, set the new members to 0 */
  x->x_tables[x->x_numtabs].x_vec = 0;
  x->x_tables[x->x_numtabs].x_arrayname = s;
  x->x_tables[x->x_numtabs].x_npoints = 0;
  x->x_whichtab = x->x_numtabs;
  x->x_numtabs++;
}

static void tabreceive_mult_remove(t_tabreceive_mult *x, t_symbol *s)
{
  int i;

  // check if the array already exists in the list
  for(i = 0; i < x->x_numtabs; i++) {
    if (x->x_tables[i].x_arrayname == s) {
      /* if the array name exists in the list, remove it by shifting
      the array elements by one position lower */
      memmove(x->x_tables + i, x->x_tables + i + 1,
        sizeof(t_tabledata) * (x->x_numtabs - i - 1));
      x->x_numtabs--;
      return;
    }
  }
}

static void tabreceive_mult_check(t_tabreceive_mult *x)
{
  int i;

  for(i = 0; i < x->x_numtabs; i++) {
    t_garray *a;

    t_word *vec = 0;
    int points;
    x->x_tables[i].x_vec = 0;
    x->x_tables[i].x_npoints = 0;

    if (!(a = (t_garray *)pd_findbyclass(x->x_tables[i].x_arrayname, garray_class)))
    {
      pd_error(x, "tabreceive_mult~: %s: no such array",
        x->x_tables[i].x_arrayname->s_name);
      /* if the wrong array name is not removed, [tabreceive_mult~] will complain
      at the next "set" message, even if the newly set array exists
      the x_init int protects from deleting an array when a patch with both
      [tabreceive_mult~] and an array it is set to listen to is loaded
      since [tabreceive_mult~] might have been created before the array itself */
      if (x->x_init) tabreceive_mult_remove(x, x->x_tables[i].x_arrayname);
      else x->x_init = 1;
    }
    /* this is the template for the function below
    garray_getfloatwords(struct _garray *x, int *size, t_word **vec)
    it assigns the array to a, its number of points to points
    and its vector to vec */
    else if (!garray_getfloatwords(a, &points, &vec))
    {
      pd_error(x, "%s: bad template for tabreceive_mult~",
        x->x_tables[i].x_arrayname->s_name);
    }
    else {
      garray_usedindsp(a);
      x->x_tables[i].x_vec = vec;
      x->x_tables[i].x_npoints = points;
    }
  }
}

static int tabreceive_mult_settable(t_tabreceive_mult *x, t_symbol *s,
                                    int argc, t_atom *argv)
{
  int i;
  int error_message = 0;
  t_sample target = 0.0f;
  t_float time_ms = 0.0f;

  switch(argc) {
  case 0:
    break;
  case 1:
    /* we're checking for the validity of the set array after we set in the structure
    and after we set the ramp values, cause if it doesn't exist, it is removed from
    the structure. this way the ramp data will also be removed */
    tabreceive_mult_addtable(x, atom_gensym(argv));
    tabreceive_mult_setramp(x, target, time_ms);
    tabreceive_mult_check(x);
    break;
  case 2:
    for (i = 0; i < argc; i++) {
      if (!i) {
        tabreceive_mult_addtable(x, atom_gensym(argv));
      }
      else {
        tabreceive_mult_setramp(x, (t_sample)atom_getfloat(argv), time_ms);
        tabreceive_mult_check(x);
      }
      argv++;
    }
    break;
  case 3:
    for (i = 0; i < argc; i++) {
      if (!i) {
        tabreceive_mult_addtable(x, atom_gensym(argv));
      }
      else {
        if (i == 1) target = (t_sample)atom_getfloat(argv);
        else {
          tabreceive_mult_setramp(x, target, atom_getfloat(argv));
          tabreceive_mult_check(x);
        }
      }
      argv++;
    }
    break;
  default:
    error_message = 1;
  }
  return error_message;
}

static void tabreceive_mult_set (t_tabreceive_mult *x, t_symbol *s, int argc, t_atom *argv)
{
  int message_error = tabreceive_mult_settable(x, s, argc, argv);

  if (message_error) pd_error(x, "tabreceive_mult~: too many elements in message");
}

static void tabreceive_mult_dsp(t_tabreceive_mult *x, t_signal **sp)
{
  /* we need to check each array at every block cycle
  this way we can add arrays that don't yet exist and they will be added
  when they are actually created */
  tabreceive_mult_check(x);
  x->x_dspticktomsec = sp[0]->s_sr / (1000 * sp[0]->s_n);
  x->x_1overn = 1./sp[0]->s_n;
  dsp_add(tabreceive_mult_perform, 3, x, sp[0]->s_vec, sp[0]->s_n);
}

static void tabreceive_mult_free(t_tabreceive_mult *x) {
  /* this is called by Pd when the object is deleted */
  freebytes(x->x_tables,
    sizeof(t_tabledata) * x->x_numtabs);
}
static void *tabreceive_mult_new(t_symbol *s, int argc, t_atom *argv)
{
  t_tabreceive_mult *x = (t_tabreceive_mult *)pd_new(tabreceive_mult_class);
  /* initialize the data structure and number of tables to 0 */
  x->x_tables = 0;
  x->x_numtabs = 0;
  x->x_init = 0;
  int argument_error = tabreceive_mult_settable(x, s, argc, argv);

  if (argument_error) pd_error(x, "tabreceive_mult~: too many arguments");

  outlet_new(&x->x_obj, &s_signal);
  return (x);
}

void tabreceive_mult_tilde_setup(void)
{
  tabreceive_mult_class = class_new(gensym("tabreceive_mult~"),
  (t_newmethod)tabreceive_mult_new, (t_method)tabreceive_mult_free,
  sizeof(t_tabreceive_mult), 0, A_GIMME, 0);
  class_addmethod(tabreceive_mult_class, (t_method)tabreceive_mult_dsp,
    gensym("dsp"), A_CANT, 0);
  class_addmethod(tabreceive_mult_class, (t_method)tabreceive_mult_set,
    gensym("set"), A_GIMME, 0);
  class_addmethod(tabreceive_mult_class, (t_method)tabreceive_mult_remove,
    gensym("remove"), A_SYMBOL, 0);
}
