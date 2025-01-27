PRO SUBCUBES, field, band, chip_nr

 data_path =  '/home/data/GNS/2015/'+band+'/' + strn(field) + '/aligned/' 
 out_path = '/home/data/GNS/2015/'+band+'/' + strn(field) + '/subcubes/'
 mask_path = data_path
 x_cube = 600 ; xaxis length of sub-cube
 y_cube = 600 ; yacis length of sub-cube
 x_large = 2400  ; xaxis length of large cube
 y_large = 1200  ; yaxis length of large cube
 valid_frac = 0.3 ; minimum fraciton of valid pixels required
 
 x_sub_shift = x_cube/2
 y_sub_shift = y_cube/2
 nx = x_large/x_sub_shift - 1
 ny = y_large/y_sub_shift - 1
 n_grid = nx * ny
 npix_sub = (long(x_cube) * long(y_cube))
;stop
 for chip = chip_nr, chip_nr do begin

  chip_nr = strn(chip)

  ; First, remove all *fits and *txt  files in target directory
  spawn, 'rm -f ' +  out_path +  'chip' + chip_nr + '/*.fits.gz'
  spawn, 'rm -f ' +  out_path +  'chip' + chip_nr + '/*.fits.gz'
  spawn, 'rm -f ' +  out_path +  'chip' + chip_nr + '/*.txt'

  readcol, data_path + 'list_chip' + chip_nr + '.txt', cube_names, FORMAT='A'
  readcol, data_path + 'masklist_chip' + chip_nr + '.txt', mask_names, FORMAT='A'
  n_cubes = n_elements(cube_names)

  for ic = 0, n_cubes - 1 do begin

   icnum = strn(ic+1)
   masks = readfits(data_path + mask_names[ic])
   cube = readfits(data_path + cube_names[ic])

   sz = size(cube)
   n3 = sz[3]

   for i_x = 0, nx -1 do begin
    for i_y = 0, ny -1 do begin
  
     xlo = i_x * x_sub_shift
     xhi = xlo + x_cube - 1
     ylo = i_y * y_sub_shift
     yhi = ylo + y_cube - 1

     ; Field valid or completely/largely masked?
     out_masks = masks[xlo:xhi,ylo:yhi,*]
     wt = total(out_masks,3)
     good = where(wt eq n3,count)
     if (float(count)/npix_sub ge valid_frac) then begin

      outnam = '_' + strn(i_x) + '_' + strn(i_y)
      out_cube = cube[xlo:xhi,ylo:yhi,*]
      writefits, out_path +  'chip' + chip_nr + '/cube' + outnam + '_' + icnum + '.fits.gz', out_cube, /COMPRESS
      writefits, out_path +  'chip' + chip_nr + '/masks'+ outnam + '_' + icnum + '.fits.gz', out_masks, /COMPRESS      
      openw, outc, out_path + 'chip' + chip_nr + '/list' + outnam + '.txt', /get_lun, /APPEND
      openw, outm, out_path  + 'chip' + chip_nr + '/masklist' + outnam + '.txt', /get_lun, /APPEND

      printf, outc,  'cube' + outnam + '_' + icnum + '.fits.gz'
      printf, outm,  'masks' + outnam + '_' + icnum + '.fits.gz'
      free_lun, outc, outm

     endif

     print, 'Done cube ' + icnum + ', sub-field ' + strn(i_x) + ' ' + strn(i_y)

    endfor
   endfor

  endfor

endfor

END
