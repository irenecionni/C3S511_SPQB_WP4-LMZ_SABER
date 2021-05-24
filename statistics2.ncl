; #############################################################################
undef(",")
function time_operations(field:numeric,
                         y1[1]:integer,
                         y2[1]:integer,
                         oper[1]:string,
                         opt[1]:string,
                         l_wgt[1]:logical)
;
; Arguments
;    field: a numeric array of rank 1 to 4, first dimension must be time.
;    y1: start year of the time period to be averaged (-1 for full range).
;    y2: end year of the time period to be averaged (-1 for full range).
;    oper: type of operations:
;            "extract": no average, just extract selected period.
;            "average": average.
;            "stddev": (unbiased) standard deviation.
;    opt: operation options (has no effect is oper = extract):
;           "annualclim": annual climatology.
;           "seasonalclim": seasonal climatology for the standard seasons DJF,
;                           MAM, JJA, SON.
;           "monthlyclim": monthly climatology jan-dec.
;           "yearly": time average over every year in [y1:y2].
;           [month strings]: climatology of selected (consecutive) months
;                            (e.g., "MAM", "SONDJ").
;           [1, 12]: climatology of the selected month ("1"=Jan, "2"=Feb, ...,
;                    "12"=Dec).
;    l_wgt: if True, calculate weighted average, with days-per-month as
;           weights (has no effect is opt = "extract").
;
; Return value
;    An array of the same rank as field or of rank-1, depending on oper/opt.
;
; Description
;    Performs differnt types of time average, standard deviation or extraction
;    of a selected time period. Weighted average (with days-per-month as
;    weights) can be optionally applied.
;
; Caveats
;    The weighted standard deviation is not yet implmented for all cases
;    The weighted standard deviation is calculated using the unbiased estimator
;    This should take into account missing values and exclude the w_i for
;    which the field contains only missing values. This feature is not
;    implemented yet.
;
; References
;
; Modification history
;    20190503-righi_mattia: removed obsolete option "mymm" (used only in
;                           reformat_obs, now outdated).
;    20140703-gottschaldt_klaus-dirk: added option "mymm".
;    20140312-righi_mattia: extended with standard deviation.
;    20140109-righi_mattia: written.
;
local funcname, scriptname, monthstr, date, year, month, idx1, idx2, loc_y1, \
  loc_y2, rank, subfield, weights, idx, idx_win, idx_spr, idx_sum, idx_aut, \
  mm, idx_1st, idx_arr, p1, p2, d2, arg, v1, v2, ym_in, years, nyear, ym, \
  dims, timec, FillValue, index
begin

  funcname = "time_operations"
  scriptname = "diag_scripts/shared/statistics.ncl"

  ;enter_msg(scriptname, funcname)

  ; Check arguments
  if (all(oper.ne.(/"extract", "average", "stddev"/))) then
    ;error_msg("f", scriptname, funcname, "unrecognized operation " + oper)
  end if

  ; Check for time dimension
  if (field!0.ne."time") then
    ;error_msg("f", scriptname, funcname, "the first dimension " + \
    ;          "of input is not time")
  end if

  ; Check for calendar attribute
  if (.not.isatt(field&time, "calendar")) then
    ;error_msg("f", scriptname, funcname, "time dimension of " + \
    ;          "input must have a calendar attribute")
  end if

  ; Define months string
  monthstr = "JFMAMJJASOND"
  monthstr = monthstr + monthstr

  ; Define flags
  l_ext = oper.eq."extract"
  l_avg = oper.eq."average"
  l_std = oper.eq."stddev"

  ; Calculate date from time coordinate
  date := cd_calendar(field&time, 0)
  year := date(:, 0)
  month := date(:, 1)

  ; Determine indexes for the requested time range
  if (y1.eq.-1) then
    idx1 = 0
    loc_y1 = toint(min(date(:, 0)))
  else
    idx1 = min(ind(year.eq.y1))
    loc_y1 = y1
  end if
  if (y2.eq.-1) then
    idx2 = dimsizes(field&time) - 1
    loc_y2 = toint(max(date(:, 0)))
  else
    idx2 = max(ind(year.eq.y2))
    loc_y2 = y2
  end if
  if (ismissing(idx1).or.ismissing(idx2)) then
    ;error_msg("f", scriptname, funcname, "the selected time " + \
     ;         "period is out of range")
  end if
  delete(date)
  delete(year)
  delete(month)

  ; Extract requested time range
  rank = dimsizes(dimsizes(field))
  if (rank.eq.4) then
    subfield = field(idx1:idx2, :, :, :)
  end if
  if (rank.eq.3) then
    subfield = field(idx1:idx2, :, :)
  end if
  if (rank.eq.2) then
    subfield = field(idx1:idx2, :)
  end if
  if (rank.eq.1) then
    subfield = field(idx1:idx2)
  end if

  ; Re-calculate date for subfield
  date := cd_calendar(subfield&time, 0)
  year := date(:, 0)
  month := date(:, 1)
  rank := dimsizes(dimsizes(subfield))

  ; Define weights as days-per-month
  if (l_wgt) then
    weights = days_in_month(toint(year), toint(month))
  else
    weights = tofloat(subfield&time)
    weights = 1.
  end if

  ; Extract only
  if (l_ext .and. opt.eq."") then
    ;leave_msg(scriptname, funcname)
    return(subfield)
  end if

  ; Calculate time average/standard deviation according to the opt argument

  ; Multi-year average
  if (opt.eq."annualclim") then
    if (l_avg) then
      out = dim_avg_wgt_n_Wrap(subfield, weights, 1, 0)
    end if
    if (l_std) then
      ;error_msg("f", scriptname, funcname, "feature not yet implemented")
    end if
    ;leave_msg(scriptname, funcname)
    return(out)
  end if

  ; Year average
  if (opt.eq."yearly") then
    ny = loc_y2 - loc_y1 + 1
    if (rank.eq.4) then
      out = subfield(0:ny - 1, :, :, :)  ; Copy metadata
      do yy = loc_y1, loc_y2
        idx = ind(year.eq.yy)
        if (l_avg) then
          out(yy - loc_y1, :, :, :) = \
            dim_avg_wgt_n_Wrap(subfield(idx, :, :, :), weights(idx), 1, 0)
        end if
        if (l_std) then
          ;error_msg("f", scriptname, funcname, "feature not yet implemented")
        end if
        delete(idx)
      end do
    end if
    if (rank.eq.3) then
      out = subfield(0:ny - 1, :, :)  ; Copy metadata
      do yy = loc_y1, loc_y2
        idx = ind(year.eq.yy)
        if (l_avg) then
          out(yy - loc_y1, :, :) = \
            dim_avg_wgt_n_Wrap(subfield(idx, :, :), weights(idx), 1, 0)
        end if
        if (l_std) then
          ;error_msg("f", scriptname, funcname, "feature not yet implemented")
        end if
        delete(idx)
      end do
    end if
    if (rank.eq.2) then
      out = subfield(0:ny - 1, :)  ; Copy metadata
      do yy = loc_y1, loc_y2
        idx = ind(year.eq.yy)
        if (l_avg) then
          out(yy - loc_y1, :) = \
            dim_avg_wgt_n_Wrap(subfield(idx, :), weights(idx), 1, 0)
        end if
        if (l_std) then
          ;error_msg("f", scriptname, funcname, "feature not yet implemented")
        end if
        delete(idx)
      end do
    end if
    if (rank.eq.1) then
      out = subfield(0:ny - 1)  ; Copy metadata
      do yy = loc_y1, loc_y2
        idx = ind(year.eq.yy)
        if (l_avg) then
          out(yy - loc_y1) = \
            dim_avg_wgt_Wrap(subfield(idx), weights(idx), 1)
        end if
        if (l_std) then
          ;error_msg("f", scriptname, funcname, "feature not yet implemented")
        end if
        delete(idx)
      end do
    end if
    out!0 = "year"
    delete(out&year)
    out&year = ispan(loc_y1, loc_y2, 1)
    ;leave_msg(scriptname, funcname)
    return(out)
  end if

  ; Season average
  if (opt.eq."seasonalclim") then
    idx_win = ind(month.eq.1.or.month.eq.2.or.month.eq.12)
    idx_spr = ind(month.eq.3.or.month.eq.4.or.month.eq.5)
    idx_sum = ind(month.eq.6.or.month.eq.7.or.month.eq.8)
    idx_aut = ind(month.eq.9.or.month.eq.10.or.month.eq.11)
    if (rank.eq.4) then
      out = subfield(0:3, :, :, :)
      if (l_avg) then
        out(0, :, :, :) = \
          dim_avg_wgt_n_Wrap(subfield(idx_win, :, :, :), \
                             weights(idx_win), 1, 0)
        out(1, :, :, :) = \
          dim_avg_wgt_n_Wrap(subfield(idx_spr, :, :, :), \
                             weights(idx_spr), 1, 0)
        out(2, :, :, :) = \
          dim_avg_wgt_n_Wrap(subfield(idx_sum, :, :, :), \
                             weights(idx_sum), 1, 0)
        out(3, :, :, :) = \
          dim_avg_wgt_n_Wrap(subfield(idx_aut, :, :, :), \
                             weights(idx_aut), 1, 0)
      end if
      if (l_std) then
        ;error_msg("f", scriptname, funcname, "feature not yet implemented")
      end if
    end if
    if (rank.eq.3) then
      out = subfield(0:3, :, :)
      if (l_avg) then
        out(0, :, :) = \
          dim_avg_wgt_n_Wrap(subfield(idx_win, :, :), weights(idx_win), 1, 0)
        out(1, :, :) = \
          dim_avg_wgt_n_Wrap(subfield(idx_spr, :, :), weights(idx_spr), 1, 0)
        out(2, :, :) = \
          dim_avg_wgt_n_Wrap(subfield(idx_sum, :, :), weights(idx_sum), 1, 0)
        out(3, :, :) = \
          dim_avg_wgt_n_Wrap(subfield(idx_aut, :, :), weights(idx_aut), 1, 0)
      end if
      if (l_std) then
        ;error_msg("f", scriptname, funcname, "feature not yet implemented")
      end if
    end if
    if (rank.eq.2) then
      out = subfield(0:3, :)
      if (l_avg) then
        out(0, :) = \
          dim_avg_wgt_n_Wrap(subfield(idx_win, :), weights(idx_win), 1, 0)
        out(1, :) = \
          dim_avg_wgt_n_Wrap(subfield(idx_spr, :), weights(idx_spr), 1, 0)
        out(2, :) = \
          dim_avg_wgt_n_Wrap(subfield(idx_sum, :), weights(idx_sum), 1, 0)
        out(3, :) = \
          dim_avg_wgt_n_Wrap(subfield(idx_aut, :), weights(idx_aut), 1, 0)
      end if
      if (l_std) then
        ;error_msg("f", scriptname, funcname, "feature not yet implemented")
      end if
    end if
    if (rank.eq.1) then
      out = subfield(0:3)
      if (l_avg)  then
        out(0) = dim_avg_wgt_Wrap(subfield(idx_win), weights(idx_win), 1)
        out(1) = dim_avg_wgt_Wrap(subfield(idx_spr), weights(idx_spr), 1)
        out(2) = dim_avg_wgt_Wrap(subfield(idx_sum), weights(idx_sum), 1)
        out(3) = dim_avg_wgt_Wrap(subfield(idx_aut), weights(idx_aut), 1)
      end if
      if (l_std) then
        ;error_msg("f", scriptname, funcname, "feature not yet implemented")
      end if
    end if
    delete(out&time)
    out!0 = "season"
    out&season = (/"DJF", "MAM", "JJA", "SON"/)
    ;leave_msg(scriptname, funcname)
    return(out)
  end if

  ; Annual cycle
  if (opt.eq."monthlyclim") then
    if (rank.eq.4) then
      out = subfield(0:11, :, :, :)  ; Copy metadata
      do mm = 0, 11
        if (l_avg)  then
          out(mm, :, :, :) = \
            dim_avg_wgt_n_Wrap(subfield(mm::12, :, :, :), \
                               weights(mm::12), 1, 0)
        end if
        if (l_std) then
          ;error_msg("f", scriptname, funcname, "feature not yet implemented")
        end if
      end do
    end if
    if (rank.eq.3) then
      out = subfield(0:11, :, :)  ; Copy metadata
      do mm = 0, 11
        if (l_avg)  then
          out(mm, :, :) = \
            dim_avg_wgt_n_Wrap(subfield(mm::12, :, :), weights(mm::12), 1, 0)
        end if
        if (l_std) then
          ;error_msg("f", scriptname, funcname, "feature not yet implemented")
        end if
      end do
    end if
    if (rank.eq.2) then
      out = subfield(0:11, :)  ; Copy metadata
      do mm = 0, 11
        if (l_avg) then
          out(mm, :) = \
            dim_avg_wgt_n_Wrap(subfield(mm::12, :), weights(mm::12), 1, 0)
        end if
        if (l_std) then
          ;error_msg("f", scriptname, funcname, "feature not yet implemented")
        end if
      end do
    end if
    if (rank.eq.1) then
      out = subfield(0:11)  ; Copy metadata
      do mm = 0, 11
        if (l_avg) then
          out(mm) = dim_avg_wgt_Wrap(subfield(mm::12), weights(mm::12), 1)
        end if
        if (l_std) then
          out(mm) = dim_stddev_wgt_Wrap(subfield(mm::12), weights(mm::12), 1)
        end if
      end do
    end if
    out!0 = "month"
    delete(out&month)
    out&month = (/"J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"/)

    ;leave_msg(scriptname, funcname)
    return(out)
  end if

  ; Months string (at least 2 consecutive months): define indexes
  if (.not.ismissing(str_match_ind_ic(monthstr, opt)).and. \
      strlen(opt).ge.2.and.strlen(opt).le.12) then
    idx_1st =  str_index_of_substr(monthstr, str_upper(opt), 1)
    idx_arr = new(strlen(opt), integer)
    do ii = 0, strlen(opt) - 1
      idx_arr(ii) = idx_1st + ii
    end do
    idx_arr = where(idx_arr.ge.12, idx_arr - 12, idx_arr)  ; Periodicity
    idx_arr = idx_arr + 1  ; From 0-based to month number
    do ii = 0, dimsizes(idx_arr) - 1
      if (.not.isdefined("idx"))  then
        idx = ind(month.eq.idx_arr(ii))
      else
        tmp = array_append_record(idx, ind(month.eq.idx_arr(ii)), 0)
        delete(idx)
        idx = tmp
        delete(tmp)
      end if
    end do
    delete(idx_1st)
    delete(idx_arr)
  end if

  ; Specific-month average: define indexes
  if (any(opt.eq.tostring(ispan(1, 12, 1)))) then
    idx = ind(month.eq.toint(opt))
  end if

  ; Extract or average over the above indexes
  if (isdefined("idx")) then
    if (rank.eq.4) then
      if (l_ext) then
        out = subfield(idx, :, :, :)
      end if
      if (l_avg) then
        if (dimsizes(idx) .eq. 1) then
          out = \
            dim_avg_wgt_n_Wrap(subfield(idx:idx, :, :, :), weights(idx), 1, 0)
        else
          out = dim_avg_wgt_n_Wrap(subfield(idx, :, :, :), weights(idx), 1, 0)
        end if
      end if
      if (l_std) then
        ;error_msg("f", scriptname, funcname, "feature not yet implemented")
      end if
    end if
    if (rank.eq.3) then
      if (l_ext) then
        out = subfield(idx, :, :)
      end if
      if (l_avg) then
        if (dimsizes(idx) .eq. 1) then
          out = dim_avg_wgt_n_Wrap(subfield(idx:idx, :, :), weights(idx), 1, 0)
        else
          out = dim_avg_wgt_n_Wrap(subfield(idx, :, :), weights(idx), 1, 0)
        end if
      end if
      if (l_std) then
        ;error_msg("f", scriptname, funcname, "feature not yet implemented")
      end if
    end if
    if (rank.eq.2) then
      if (l_ext) then
        out = subfield(idx, :)
      end if
      if (l_avg) then
        if (dimsizes(idx) .eq. 1) then
          out = dim_avg_wgt_n_Wrap(subfield(idx:idx, :), weights(idx), 1, 0)
        else
          out = dim_avg_wgt_n_Wrap(subfield(idx, :), weights(idx), 1, 0)
        end if
      end if
      if (l_std) then
        ;error_msg("f", scriptname, funcname, "feature not yet implemented")
      end if
    end if
    if (rank.eq.1) then
      if (l_ext) then
        out = subfield(idx)
      end if
      if (l_avg) then
        if (dimsizes(idx) .eq. 1) then
          out = dim_avg_wgt_n_Wrap(subfield(idx:idx), weights(idx), 1, 0)
        else
          out = dim_avg_wgt_n_Wrap(subfield(idx), weights(idx), 1, 0)
        end if
      end if
      if (l_std) then
        ;error_msg("f", scriptname, funcname, "feature not yet implemented")
      end if
    end if
    ;leave_msg(scriptname, funcname)
    return(out)
  end if

  ;error_msg("f", scriptname, funcname, "unrecognized option " + opt)

end
undef("coswgt_area_average")
function coswgt_area_average(field:numeric)
;
; Arguments
;    field: numeric field.
;
; Return value
;    The area average using cosine lat weights.
;
; Description
;    Computes the area average using cosine lat weights and lon weights=1.
;
local funcname, scriptname, lat, wgt_lat, lon, lon_size, wgt_lon, ave
begin

  funcname = "coswgt_area_average"
  scriptname = "diag_scripts/shared/statistics.ncl"

  lat = field&lat
  wgt_lat = tofloat(NormCosWgtGlobe(lat))

  lon = field&lon
  lon_size = dimsizes(lon)
  wgt_lon = new((/lon_size(0)/), float)
  wgt_lon = 1.0

  ave = wgt_areaave_Wrap(field, wgt_lat, wgt_lon, 0)
  return(ave)

end
; #############################################################################
undef("coswgt_area_std")
function coswgt_area_std(field:numeric)
;
; Arguments
;    field: numeric field.
;
; Return value
;    The area average using cosine lat weights.
;
; Description
;    Computes the area average using cosine lat weights and lon weights=1.
;
begin

  funcname = "coswgt_area_std"

   lat = field&lat
   lon = field&lon
   re   = 6.37122e06
   rad  = 4.0 * atan(1.0) / 180.0
   con  = re * rad
   clat = cos(lat * rad)           ; cosine of latitude
   dlon = (lon(2) - lon(1))        ; assume dlon is constant
   dlat = (lat(2) - lat(1))        ; assume dlat is constant

   dx   = con * dlon * clat        ; dx at each latitude
   dy   = con * dlat               ; dy is constant
   dydx = dy * dx                  ; dydx(nlat)
   nlat=dimsizes(lat)
   mlon=dimsizes(lon)
   wgt  = new((/nlat, mlon/), typeof(dydx))
   wgt=conform (wgt, dydx, 0)



   ave = wgt_areaave2(field, wgt, 0)
   d2 = (field - ave) ^ 2 ;
     ; Calculate weighted standard deviation
   ss=sum(d2)
   arg = wgt_areasum2(d2, wgt, 0)

   v1 = sum(where(ismissing(field), wgt@_FillValue, wgt))
   v2 = sum(where(ismissing(field), wgt@_FillValue, wgt) ^ 2)
   dd=(v1 / (v1 ^ 2 - v2))

   out = sqrt( arg*dd)

  return(out)


end
  ; Pi [1]
    PI = 3.1415927
    PI@units = "1"
    PI@name = "Greek Pi"


    EARTHradius = 6371000.
    EARTHradius@units = "m"
    EARTHradius@name = "Earth radius"
    ; Earth surface area [m^2]
  
    EARTHarea = 4. * PI * EARTHradius ^ 2
    EARTHarea@units = "m2"
    EARTHarea@name = "Earth surface area"
; #############################################################################
undef("gridcell_area")
function gridcell_area(deltax[1]: numeric,
                       lat_lo[1]: numeric,
                       lat_hi[1]: numeric)
local funcname, scriptname, deg2rad
begin

  funcname = "gridcell_area"
  scriptname = "diag_scripts/shared/latlon.ncl"
  ;enter_msg(scriptname, funcname)

  ; Convert to radiants
  deg2rad = acos(-1.0) / 180.

  deltay = sin(lat_hi * deg2rad) - sin(lat_lo * deg2rad)
  out = abs(EARTHradius ^ 2 * deltay * deltax * deg2rad)
  out@units = EARTHarea@units

  ;leave_msg(scriptname, funcname)
  return(out)

end

; #############################################################################
; #############################################################################
undef("map_area")
function map_area(lat[*]:numeric,
                  lon[*]:numeric)
local funcname, scriptname, loclat, loclon, nlat, lati, xresol
begin

  funcname = "map_area"
  scriptname = "diag_scripts/shared/latlon.ncl"
  ;enter_msg(scriptname, funcname)

  ; Type conversion
  loclat = tofloat(lat)
  loclon = tofloat(lon)
  nlat = dimsizes(loclat)

  ; Check ordering (must be S-N)
  if (isMonotonic(loclat).ne.1) then
    ;error_msg("f", scriptname, funcname, "input latitude is not " + \
;              "monotonically increasing")
  end if

  ; Calculate longitude resolution (assumed constant)
  xresol = abs(loclon(1) - loclon(0))
  dsize = dimsizes(loclon)
  if (dsize(0).gt.2) then
    if ((xresol - abs(loclon(2) - loclon(1))).gt.1.0e-4) then
      ;error_msg("f", scriptname, funcname, "input longitude " + \
 ;               "resolution is not constant")
    end if
  end if

  ; Calculate latitude interfaces
  lati = new(dimsizes(loclat) + 1, float)
  lati(0) = max((/(3 * loclat(0) - loclat(1)) / 2., -90./))
  do ii = 1, dimsizes(lati) - 2
    lati(ii) = 0.5 * (loclat(ii - 1) + loclat(ii))
  end do
  lati(dimsizes(lati) - 1) = \
    min((/(3 * loclat(nlat - 1) - loclat(nlat - 2)) / 2., 90./))

  ; Calculate area
  area = new((/dimsizes(loclat), dimsizes(loclon)/), float)
  do ii = 0, dimsizes(loclat) - 1
    area(ii, :) = gridcell_area(xresol, lati(ii), lati(ii + 1))
  end do
  ; Assign coordinates
  area!0   = "lat"
  area&lat = lat
  area!1   = "lon"
  area&lon = lon

  ;leave_msg(scriptname, funcname)
  return(area)



end
;************************************************
function area_operations(field:numeric,
                         latmin[1]:numeric,
                         latmax[1]:numeric,
                         lonmin[1]:numeric,
                         lonmax[1]:numeric,
                         opt[1]:string,
                         l_wgt[1]:logical)

local funcname, scriptname, verbosity, rank, dimnames, wgt, conwgt, f_latmin, \
    f_latmax, f_lonmin, f_lonmax, delta, idx_latmin, idx_latmax, latidx, \
    idx_lonmin, idx_lonmax, lonidx, subfield, subwgt, dnames, didx_r, didx_l, \
    norm
begin

    funcname = "area_operations"
    scriptname = "diag_scripts/lib/ncl/latlon.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    ;enter_msg(scriptname, funcname, 6)

    ;; Check rank
    if (dimsizes(dimsizes(field)) .lt. 2) then
       ; error_msg("f", scriptname, funcname, "rank of input field <2")
    end if

    ;; Check dimensions
    rank = dimsizes(dimsizes(field))
  dimnames = getvardims(field)
    if (dimnames(rank - 1) .ne. "lon") then
        ;error_msg("f", scriptname, funcname, "the rightmost dimension " + \
              ;    "of the input field must be lon")
    end if
    if (dimnames(rank - 2) .ne. "lat") then
        ;error_msg("f", scriptname, funcname, "the second rightmost " + \
               ;   "dimension of the input field must be lat")
    end if
    delete(dimnames)

    ;; Check that grid is regular
    if (dimsizes(dimsizes(field&lat)) .ne. 1 .or. \
        dimsizes(dimsizes(field&lon)) .ne. 1) then
        ;error_msg("f", scriptname, funcname, "irregular grids not supported")
    end if

    ;; Calculate weights
    wgt = map_area(field&lat, field&lon)
    if (.not.l_wgt) then
        wgt = 1.
    end if

    ;; Conform weights and filter for missing values in the input field
    if (rank.eq.4) then
        conwgt = conform(field, wgt, (/2, 3/))
    end if
    if (rank.eq.3) then
        conwgt = conform(field, wgt, (/1, 2/))
    end if
    if (rank.eq.2) then
        conwgt = wgt
    end if
    conwgt = where(ismissing(field), conwgt@_FillValue, conwgt)
    copy_VarCoords(field, conwgt)
    delete(wgt)

    ;; Copy arguments to local (swap lat limits if necessary)
    f_latmin = todouble(min((/latmin, latmax/)))
    f_latmax = todouble(max((/latmin, latmax/)))
    f_lonmin = todouble(where(lonmin .lt. 0., lonmin + 360., lonmin))
    f_lonmax = todouble(where(lonmax .lt. 0., lonmax + 360., lonmax))
; If single gridbox lon = 360, replace with 0 to account for periodicity
    if (f_lonmin .eq. 360 .and. f_lonmax .eq. 360) then
        f_lonmin = 0.d0
        f_lonmax = 0.d0
   end if

    ;; Find nearest points to latitude boundaries
    ;; In the (unlikely) case that the selected boundary point is at the
    ;; midpoint between two gridboxes, ind will return two points.
    ;; The minimum index (= minumum value) is chosen in such case,
    ;; consistently with NCL coordinate subscripting
    ;; http://www.ncl.ucar.edu/Document/Language/subscript_coord.shtml
    delta = abs(field&lat - f_latmin)
    idx_latmin = min(ind(delta .eq. min(delta)))
    delta = abs(field&lat - f_latmax)
    idx_latmax = min(ind(delta .eq. min(delta)))
    delete(delta)
    if (ismissing(idx_latmin) .or. ismissing(idx_latmax)) then
        ;error_msg("f", scriptname, funcname, "no values found in the " + \
               ;   "given latitude range")
 end if
    if (latmin .eq. latmax .and. idx_latmin .ne. idx_latmax) then
        ;error_msg("f", scriptname, funcname, "multiple indexes found " + \
               ;   "for single gridbox latitude selection")
    end if

    ;; Define latitude indexes
    latidx = ispan(idx_latmin, idx_latmax, 1)  ; array of indexes
    delete(idx_latmin)
    delete(idx_latmax)

    ;; Find neartest points to longitude boundaries
    ;; In the (unlikely) case that the selected boundary point is at the
    ;; midpoint between two gridboxes, ind will return two points.
    ;; The minimum index (= minumum value) is chosen in such case,
    ;; consistently with NCL coordinate subscripting
    ;; http://www.ncl.ucar.edu/Document/Language/subscript_coord.shtml
    if (f_lonmin.gt.max(field&lon)) then
        idx_lonmin = 0
    else
        delta = abs(field&lon - f_lonmin)
        idx_lonmin = min(ind(delta .eq. min(delta)))
        delete(delta)
    end if
    if (f_lonmax.gt.max(field&lon)) then
        idx_lonmax = dimsizes(field&lon) - 1
    else
        delta = abs(field&lon - f_lonmax)
        idx_lonmax = min(ind(delta .eq. min(delta)))
        delete(delta)
    end if
    if (ismissing(idx_lonmin) .or. ismissing(idx_lonmax)) then
        ;error_msg("f", scriptname, funcname, "no values found in the " + \
                  ;"given longitude range")
    end if
    ;; Special case: single gridbox selection beteen max(field&lon) and 0
    if (lonmin .eq. lonmax .and. lonmin.gt.max(field&lon)) then
        idx_lonmin = where(abs(360. - lonmin) .lt. \
                           abs(max(field&lon) - lonmin), \
                           0, dimsizes(field&lon) - 1)
        idx_lonmax = idx_lonmin
    end if
    if (lonmin .eq. lonmax .and. idx_lonmin .ne. idx_lonmax) then
        ;error_msg("f", scriptname, funcname, "multiple indexes found " + \
                  ;"for single gridbox longitude selection")
    end if

    ;; Define longitude indexes
    if (idx_lonmin .le. idx_lonmax) then  ; Greenwich outside
        lonidx = ispan(idx_lonmin, idx_lonmax, 1)
    else  ; Greenwich inside
        lonidx = array_append_record( \
            ispan(idx_lonmin, dimsizes(field&lon) - 1, 1), \
            ispan(0, idx_lonmax, 1), 0)
    end if
    delete(idx_lonmin)
    delete(idx_lonmax)
;; Extract region
    if (rank .eq. 4) then
        subfield = field(:, :, latidx, lonidx)
        subwgt = conwgt(:, :, latidx, lonidx)
    end if
    if (rank .eq. 3) then
        subfield = field(:, latidx, lonidx)
        subwgt = conwgt(:, latidx, lonidx)
    end if
    if (rank .eq. 2) then
        subfield = field(latidx, lonidx)
      subwgt = conwgt(latidx, lonidx)
    end if
    delete(conwgt)

    if (opt .eq. "extract") then
        ;;leave _msg(scriptname, funcname, 6)
        return(subfield)
    end if

    ;; Find dimensions to be averaged or summed
    dnames = getVarDimNames(subfield)
    didx_r = ind(dnames.eq."lat" .or. dnames.eq."lon")
    didx_l = ind(dnames.ne."lat" .and. dnames.ne."lon")
    if (all(ismissing(didx_r))) then
        ;error_msg("f", scriptname, funcname, "cannot perform " + \
        ;          "average/sum on single gridcell, use the extract option")
    end if

    ;; Normalize weights
    if (opt .eq. "average") then
        if (dimsizes(didx_r).eq.dimsizes(dnames)) then
            subwgt = subwgt / sum(subwgt)
        else
            norm = dim_sum_n_Wrap(subwgt, didx_r)
            subwgt = subwgt / conform(subwgt, norm, didx_l)
        end if
    end if

    ;; Calculate average/sum. The same code applies, as the only difference
    ;; is in the normalization of the weights, which is calculated above
    ;; only for opt="average"
    if (opt .eq. "average" .or. opt .eq. "sum") then
        subfield = subfield * subwgt
        if (dimsizes(didx_r).eq.dimsizes(dnames)) then
            out = dim_sum_Wrap(dim_sum_Wrap(subfield))
        else
            out = dim_sum_n_Wrap(subfield, didx_r)
        end if
       ; ;leave _msg(scriptname, funcname, 6)
        return(out)
    end if
  ;error_msg("f", scriptname, funcname, "unrecognized option " + opt)

end

