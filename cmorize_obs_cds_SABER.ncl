; #############################################################################
;  CMORizer for SABER data
; #############################################################################
; ############################################################################
load "/Users/irenecionni/Copernicus/utilities.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
begin

  ; Script name (for logger)
  DIAG_SCRIPT = "cmorize_obs_cds_saber.ncl"

  ; Source name
  OBSNAME = "SABER"

  ; Tier
  TIER = 3

  ; Period
  YEAR1 = 2002
  YEAR2 = 2020

  ; Selected variable (standard name)
  VAR = (/"tro3"/)

  ; Name in the raw data
  ;NAME = (/"ozone_mixing_ratio"/)

  ; MIP
  MIP = (/"Amon"/)

  ; Frequency
  FREQ = (/"mon"/)

  ; CMOR table
  CMOR_TABLE = "CMOR_tro3.dat" 
  ; Type
  TYPE = "sat"

  ; Version
  VERSION = "L3"

  ; Global attributes
  SOURCE = "https://cds.climate.copernicus.eu/cdsapp#!/dataset/satellite-ozone?tab=form"
  REF = ""
  COMMENT = ""

end

begin

   input_dir_path="/dataset-satellite-ozone-SABER/"
   output_dir_path="/SABER/"
  ; Loop over variables
  do vv = 0, dimsizes(VAR) - 1
    nlat=18
    nlev=26
    nlon=36

    yyyymm = yyyymm_time(YEAR1,YEAR2,"integer")
    print(yyyymm)
    yyyy  = yyyymm/100 
  mm  = yyyymm-yyyy*100                ; mmdd = yyyymmdd%10000
  dd    = mm                        ; dd   = mmdd%100

  hh    = mm                                 ; create arrays [*] of required size
  mn    = mm
  sc    = mm
  dd    = 15
  hh    = 0                                  ; array syntax
  mn    = 0
  sc    = 0
   units  = "hours since 1900-01-01 00:00:00"
    
 time=cd_inv_calendar(yyyy,mm,dd,hh,mn,sc,units, 0)
 time!0 = "time"
  print(time)

    output1=new((/dimsizes(yyyymm),nlev,nlat,nlon/),"double")
    do iy=0, dimsizes(yyyymm)-1
       print(yyyymm(iy))

       fname = systemfunc("ls "+ input_dir_path + yyyymm(iy)+"-C3S-L3_OZONE-O3_PRODUCTS-SABER-TIMED-ALG-MONTHLY-v0002.nc")
       print(fname)
       ;setfileoption("nc", "MissingToFillValue", False)
       if (.not.ismissing(fname)) then
          f = addfile(fname, "r")
       ;ListSetType (f, "join")
       do il=0,nlon-1 
         output1(iy,:,:,il) =(/ f->ozone_mixing_ratio/)
       end do
         lat=f->latitude_centers
         plev=f->pressure
       end if
    end do
 ;   time=time_or(: dimsizes(yyyymm)-1)
    print(time)
    print(yyyymm(: dimsizes(yyyymm)-1))
    printVarSummary(output1)
    output1!0="time"
    output1&time=time
    output1!1="plev"
    output1&plev=plev
     output1!2="lat" 
    output1&lat=lat
         output1!3="lon"
   lon=lonGlobeF(nlon, "lon", "longitude", "degrees_east")
   output1&lon=lon
 
    printVarSummary(output1)
    output1!0="time"
     output=output1
    output=(/output1/)
    printVarSummary(output) 
   ; Format coordinates
    
    format_coords(output, YEAR1 + "0101", YEAR2 + "1231", FREQ(vv))
     printVarSummary(output)
    ; Set variable attributes
    tmp = format_variable(output, VAR(vv), CMOR_TABLE(vv))
    delete(output)
    output = tmp
    printVarSummary(tmp)
    print(tmp&time)
    delete(tmp)
    ; Calculate coordinate bounds
    bounds = guess_coord_bounds(output, FREQ(vv))

    ; Set global attributes
    ;gAtt = set_global_atts(OBSNAME, TIER, SOURCE, REF, COMMENT, CMOR_TABLE(vv))
    gAtt = set_global_atts(OBSNAME, TIER, SOURCE, REF, COMMENT)
    ; Output file
    DATESTR = YEAR1 + "01-" + YEAR2 + "12"
    fout = output_path + \
      str_join((/"OBS", OBSNAME, TYPE, VERSION, \
                 MIP(vv), VAR(vv), DATESTR/), "_") + ".nc"

    ; Write variable
    print(bounds)
    write_nc(fout, VAR(vv), output, bounds, gAtt)
    delete(gAtt)
    delete(output)
    delete(bounds)

  end do

end
