
/****************************************************************************
 *
 * MODULE:       r.heliosat
 * AUTHOR(S):    Quinn Hart
 * PURPOSE:      Borrowed from r.sunhours, this application computes, NRELs solpos on
 *               a raster, and can provide solar azimuth, angle, sunshine hours, sunrise/sunset time,
 *               and an instantaneous sun angle.
 *               Uses NREL SOLPOS
 * COPYRIGHT:    (C) 2018 by the GRASS Development Team
 *
 *               This program is free software under the GNU General Public
 *               License (>=v2). Read the file COPYING that comes with GRASS
 *               for details.
 *
 *****************************************************************************/

#define MAIN

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <grass/gis.h>
#include <grass/raster.h>
#include <grass/gprojects.h>
#include <grass/glocale.h>
#include "solpos00.h"

//Taken from heliosat code.
// Define the value of PI
#define PI 3.141592653589793

/* internal undefined value for NULL */
#define UNDEFZ      -9999.

// Coefficients for the beam angular function, Fb
const double AP[3][3] = {
    { 2.6463e-01, -6.1581e-02,  3.1408e-03 },
    { 2.0402    ,  1.8945e-02, -1.1161e-02 },
    { -1.3025   ,  3.9231e-02,  8.5079e-03 }
};

const double LE[3][3][4] = {
    // y > 30
    {
        { -1.7349e-2, -5.8985e-3,  6.8868e-4, 0 },
        {  1.0258   , -1.2196e-1,  1.9229e-3, 0 },
        { -7.2178e-3,  1.3086e-1, -2.8405e-3, 0 }
    },
    // 15 < y < 30
    {
        { -8.2193e-3,  4.5643e-4,  6.7916e-5, 0 },
        {  8.9233e-1, -1.9991e-1,  9.9741e-3, 0 },
        {  2.5428e-1,  2.6140e-1, -1.7020e-2, 0 }
    },
    // y < 15
    {
        { -1.1656e-3,  1.8408e-4, -4.8754e-7, 0 },
        {  7.4095e-1, -2.2427e-1,  1.5314e-2, 0 },
        {  3.4959e-1,  7.2313e-1, -1.2305e-1, 5.9194e-3 }
    }
};

// Right now use only Y > 30
const double L[3][3] = {
    { -1.7349e-2, -5.8985e-3,  6.8868e-4 },
    {  1.0258   , -1.2196e-1,  1.9229e-3 },
    { -7.2178e-3,  1.3086e-1, -2.8405e-3 }
};



void set_solpos_time(struct posdata *pdat, int year, int month, int day, int hour, int minute, int second, int timezone);
void set_solpos_longitude(struct posdata *, double );
int roundoff(double *);


int main(int argc, char *argv[])
{
    struct GModule *module;
    struct {
      struct Option *elevin, *linkein, *beam, *diffuse, *total, *sunhours, *year, *month, *day, *hour, *minutes, *seconds, *timezone;
    } parm;
    struct Cell_head window;
    FCELL *elevinbuf, *linkeinbuf, *beambuf, *diffusebuf, *totalbuf, *sunhourbuf;
    struct History hist;

    /* projection information of input map */
    struct Key_Value *in_proj_info, *in_unit_info;
    struct pj_info iproj; /* input map proj parameters  */
    struct pj_info oproj; /* output map proj parameters  */
    struct pj_info tproj; /* transformation parameters  */
    char *elevin_name, *linkein_name, *beam_name, *diffuse_name, *total_name, *sunhour_name;
    int elevin_fd, linkein_fd, beam_fd, diffuse_fd, total_fd, sunhour_fd;
    double east, east_ll, north, north_ll;
    double north_gc, north_gc_sin, north_gc_cos; /* geocentric latitude */
    double ba2;
    int year, month, day, hour, minutes, seconds, timezone;
    int doy; /* day of year */
    int row, col, nrows, ncols;
    int do_reproj = 0;
    struct posdata pd;

    G_gisinit(argv[0]);

    module = G_define_module();
    G_add_keyword(_("raster"));
    G_add_keyword(_("solar"));
    G_add_keyword(_("sun energy"));
    G_add_keyword(_("sun position"));
    module->label =
        _("Calculates daily integrated radiation parameters on a raster map");
    module->description =
        _("beam: direct beam radiation"
          "diffuse: diffuse radiation"
          "total: total radiation");

    parm.elevin = G_define_option();
    parm.elevin->key = "elevation";
    parm.elevin->type = TYPE_STRING;
    parm.elevin->required = YES;
    parm.elevin->gisprompt = "old,cell,raster";
    parm.elevin->description =
        _("Name of the input elevation raster map [meters]");
    parm.elevin->guisection = _("Input");

    parm.linkein = G_define_option();
    parm.linkein->key = "linke";
    parm.linkein->type = TYPE_STRING;
    parm.linkein->required = NO;
    parm.linkein->gisprompt = "old,cell,raster";
    parm.linkein->description = _("Name of the Linke atmospheric turbidity "
                                  "coefficient input raster map [-]");
    parm.linkein->guisection = _("Input");

    parm.beam = G_define_standard_option(G_OPT_R_OUTPUT);
    parm.beam->key = "beam";
    parm.beam->label = _("Output raster map with integrated beam radiation");
    parm.beam->required = NO;

    parm.diffuse = G_define_standard_option(G_OPT_R_OUTPUT);
    parm.diffuse->key = "diffuse";
    parm.diffuse->label = _("Output raster map with integrated diffuse radiation");
    parm.diffuse->required = NO;

    parm.total = G_define_standard_option(G_OPT_R_OUTPUT);
    parm.total->key = "total";
    parm.total->label = _("Output raster map with integrated total radiation");
    parm.total->required = NO;

    parm.sunhours = G_define_standard_option(G_OPT_R_OUTPUT);
    parm.sunhours->key = "sunhour";
    parm.sunhours->label = _("Output raster map with sunshine hours");
    parm.sunhours->required = NO;

    parm.year = G_define_option();
    parm.year->key = "year";
    parm.year->type = TYPE_INTEGER;
    parm.year->required = YES;
    parm.year->description = _("Year");
    parm.year->options = "1950-2050";
    parm.year->guisection = _("Time");

    parm.month = G_define_option();
    parm.month->key = "month";
    parm.month->type = TYPE_INTEGER;
    parm.month->required = NO;
    parm.month->label = _("Month");
    parm.month->description =
        _("If not given, day is interpreted as day of the year");
    parm.month->options = "1-12";
    parm.month->guisection = _("Time");

    parm.day = G_define_option();
    parm.day->key = "day";
    parm.day->type = TYPE_INTEGER;
    parm.day->required = YES;
    parm.day->description = _("Day");
    parm.day->options = "1-366";
    parm.day->guisection = _("Time");

    parm.hour = G_define_option();
    parm.hour->key = "hour";
    parm.hour->type = TYPE_INTEGER;
    parm.hour->required = NO;
    parm.hour->description = _("Hour");
    parm.hour->options = "0-24";
    parm.hour->answer = "12";
    parm.hour->guisection = _("Time");

    parm.minutes = G_define_option();
    parm.minutes->key = "minute";
    parm.minutes->type = TYPE_INTEGER;
    parm.minutes->required = NO;
    parm.minutes->description = _("Minutes");
    parm.minutes->options = "0-60";
    parm.minutes->answer = "0";
    parm.minutes->guisection = _("Time");

    parm.seconds = G_define_option();
    parm.seconds->key = "second";
    parm.seconds->type = TYPE_INTEGER;
    parm.seconds->required = NO;
    parm.seconds->description = _("Seconds");
    parm.seconds->options = "0-60";
    parm.seconds->answer = "0";
    parm.seconds->guisection = _("Time");

    parm.timezone = G_define_option();
    parm.timezone->key = "timezone";
    parm.timezone->type = TYPE_INTEGER;
    parm.timezone->required = NO;
    parm.timezone->description = _("Timezone");
    parm.timezone->options = "-12-12";
    parm.timezone->answer = "0";
    parm.timezone->guisection = _("Time");

    if (G_parser(argc, argv))
      exit(EXIT_FAILURE);

    G_get_window(&window);

    /* require at least one output or report */
    beam_name = parm.beam->answer;
    diffuse_name = parm.diffuse->answer;
    total_name = parm.total->answer;
    sunhour_name = parm.sunhours->answer;
    if (!beam_name && !diffuse_name && !total_name && !sunhour_name)
      G_fatal_error(_("No output requested, exiting."));

    year = atoi(parm.year->answer);
    if (parm.month->answer)
        month = atoi(parm.month->answer);
    else
        month = -1;

    day = atoi(parm.day->answer);
    hour = atoi(parm.hour->answer);
    minutes = atoi(parm.minutes->answer);
    seconds = atoi(parm.seconds->answer);
    timezone = atoi(parm.timezone->answer);

    /* init variables */
    north_gc_cos = 0;
    north_gc_sin = 1;

    if ((G_projection() != PROJECTION_LL)) {
        if (window.proj == 0)
            G_fatal_error(_("Current projection is x,y (undefined)."));

        do_reproj = 1;

        /* read current projection info */
        if ((in_proj_info = G_get_projinfo()) == NULL)
            G_fatal_error(_("Cannot get projection info of current location"));

        if ((in_unit_info = G_get_projunits()) == NULL)
            G_fatal_error(_("Cannot get projection units of current location"));

        if (pj_get_kv(&iproj, in_proj_info, in_unit_info) < 0)
            G_fatal_error(
                _("Cannot get projection key values of current location"));

        G_free_key_value(in_proj_info);
        G_free_key_value(in_unit_info);

        oproj.pj = NULL;
        tproj.def = NULL;

        if (GPJ_init_transform(&iproj, &oproj, &tproj) < 0)
            G_fatal_error(_("Unable to initialize coordinate transformation"));
    }

    /* always init pd */
    S_init(&pd);

    if (month == -1)
        doy = day;
    else
        doy = dom2doy2(year, month, day);

    set_solpos_time(&pd, year, 1, doy, hour, minutes, seconds, timezone);
    set_solpos_longitude(&pd, 0);
    pd.latitude = 0;

    ba2 = 6356752.3142 / 6378137.0;
    ba2 = ba2 * ba2;

    // solar functions needed
    pd.function = S_GEOM | S_ZENETR | S_SRSS | S_ETR | S_SSHA ;

    S_solpos(&pd);

    if (beam_name) {
        if ((beam_fd = Rast_open_new(beam_name, FCELL_TYPE)) < 0)
            G_fatal_error(_("Unable to create raster map <%s>"), beam_name);

        beambuf = Rast_allocate_f_buf();
    }
    else {
        beambuf = NULL;
        beam_fd = -1;
    }

    if (diffuse_name) {
      if ((diffuse_fd = Rast_open_new(diffuse_name, FCELL_TYPE)) < 0)
        G_fatal_error(_("Unable to create raster map <%s>"), diffuse_name);

      diffusebuf = Rast_allocate_f_buf();
    }
    else {
      diffusebuf = NULL;
      diffuse_fd = -1;
    }

    if (total_name) {
      if ((total_fd = Rast_open_new(total_name, FCELL_TYPE)) < 0)
        G_fatal_error(_("Unable to create raster map <%s>"), total_name);

      totalbuf = Rast_allocate_f_buf();
    }
    else {
      totalbuf = NULL;
      total_fd = -1;
    }

    if (sunhour_name) {
        if ((sunhour_fd = Rast_open_new(sunhour_name, FCELL_TYPE)) < 0)
            G_fatal_error(_("Unable to create raster map <%s>"), sunhour_name);

        sunhourbuf = Rast_allocate_f_buf();
    }
    else {
        sunhourbuf = NULL;
        sunhour_fd = -1;
    }

    nrows = Rast_window_rows();
    ncols = Rast_window_cols();

    elevin_fd = Rast_open_old(parm.elevin->answer, "");
    if (elevin_fd < 0)
        G_fatal_error(_("Unable to open raster map <%s>"), parm.elevin->answer);

    elevinbuf = Rast_allocate_f_buf();

    linkein_fd = Rast_open_old(parm.linkein->answer, "");
    if (linkein_fd < 0)
      G_fatal_error(_("Unable to open raster map <%s>"), parm.linkein->answer);

    linkeinbuf = Rast_allocate_f_buf();


    for (row = 0; row < nrows; row++) {

        G_percent(row, nrows, 2);
        // elevation && linke turbidity
        Rast_get_f_row(elevin_fd, elevinbuf, row);
        Rast_get_f_row(linkein_fd, linkeinbuf, row);

        /* get cell center northing */
        north = window.north - (row + 0.5) * window.ns_res;
        north_ll = north;

        for (col = 0; col < ncols; col++) {
            long int retval;
            double z, linke;
            // Heliosat variables
            double A0, A1, A2, B0, B1, B2, Bc, Bci, Bcz, C0, C1, C2;
            double D0, D1, D2, Dc, Dci, Dcz, Fbi, Fbi_sr, Fbiss, Fdi, Fdi_sr, Fdiss;
            double Gc, Gi, Trb, Trd, clcd, ha, hclinke, p_p0, ray, slsd;
            bool is_not_null = true;

            if (!Rast_is_f_null_value(elevinbuf + col))
              z = (float)elevinbuf[col];
            else
              is_not_null = false;

            if (!Rast_is_f_null_value(linkeinbuf + col))
              linke = (float)linkeinbuf[col];
            else
              is_not_null = false;

            if (is_not_null) {
              /* get cell center easting */
              east = window.west + (col + 0.5) * window.ew_res;
              east_ll = east;

              if (do_reproj) {
                north_ll = north;
                if (GPJ_transform(&iproj, &oproj, &tproj, PJ_FWD, &east_ll,
                                  &north_ll, NULL) < 0)
                  G_fatal_error(
                        _("Error in %s (projection of input coordinate pair)"),
                        "GPJ_transform()");
              }

              /* geocentric latitude */
              north_gc = atan(ba2 * tan(DEG2RAD * north_ll));
              north_gc_sin = sin(north_gc);
              roundoff(&north_gc_sin);
              north_gc_cos = cos(north_gc);
              roundoff(&north_gc_cos);

              set_solpos_longitude(&pd, east_ll);
              pd.latitude = north_gc * RAD2DEG;
              retval = S_solpos(&pd);
              S_decode(retval, &pd);
              G_debug(3, "solpos hour angle: %.5f", pd.hrang);

              // Common variables
              p_p0=exp(-z/8434.5);
              hclinke=linke*p_p0;
              ray=1/(6.6296+1.7513*p_p0-0.1202*pow(p_p0,2)+0.0065*pow(p_p0,3));
              slsd=sin(pd.latitude*DEG2RAD)*sin(pd.declin*DEG2RAD);
              clcd=cos(pd.latitude*DEG2RAD)*cos(pd.declin*DEG2RAD);

              // Beam Parameters
              Trb=exp(-0.8662*hclinke*ray);
              Bcz=pd.etrn*Trb;
              C0=L[0][0]+L[0][1]*hclinke+L[0][2]*pow(hclinke,2);
              C1=L[1][0]+L[1][1]*hclinke+L[1][2]*pow(hclinke,2);
              C2=L[2][0]+L[2][1]*hclinke+L[2][2]*pow(hclinke,2)+L[2][3]*pow(hclinke,3);
              B0=C0+C1*slsd+C2*pow(slsd,2)+0.5*C2*pow(clcd,2);
              B1=C1*clcd+2*C2*slsd*clcd;
              B2=0.25*C2*pow(clcd,2);
              Fbiss=B0*pd.ssha*PI/180+B1*sin(pd.ssha*DEG2RAD)+B2*sin(2*pd.ssha*DEG2RAD);
              Bc=2*Fbiss*Bcz*(12/PI);

              // Diffuse Parameters
              Trd=-1.5834e-2+3.03543e-2*hclinke+3.797e-4*pow(hclinke,2);
              Dcz=pd.etrn*Trd;
              A0=AP[0][0]+AP[0][1]*hclinke+AP[0][2]*pow(hclinke,2);
              A1=AP[1][0]+AP[1][1]*hclinke+AP[1][2]*pow(hclinke,2);
              A2=AP[2][0]+AP[2][1]*hclinke+AP[2][2]*pow(hclinke,2);
              D0=A0+A1*slsd+A2*pow(slsd,2)+0.5*A2*pow(clcd,2);
              D1=A1*clcd+2*A2*slsd*clcd;
              D2=0.25*A2*pow(clcd,2);
              Fdiss=D0*pd.ssha*PI/180+D1*sin(pd.ssha*DEG2RAD)+D2*sin(2*pd.ssha*DEG2RAD);
              Dc=2*Fdiss*Dcz*(12/PI);

              // Hour angle - max and min = sunrise and sunset
              //ha=fmin(pd.ssha,fmax(-pd.ssha,pd.hrang));
              ha=pd.hrang;

              //              if (pd.tst < pd.sretr || pd.tst > pd.ssetr)
              //  G_warning(_("Solar time is outside of sunrise/sunset time"));

              if (pd.tst < pd.sretr) {
                Bci=0;
                Dci=0;
              } else if (pd.tst > pd.ssetr) {
                Bci=Bc;
                Dci=Dc;
              } else {
                // Integrated Beam Parameters
                Fbi=B0*ha*PI/180+B1*sin(ha*DEG2RAD)+B2*sin(2*ha*DEG2RAD);
                Bci=(12/PI)*(Fbi + Fbiss)*Bcz;

                // Integrated Diffuse Parameters
                Fdi=D0*ha*PI/180+D1*sin(ha*DEG2RAD)+D2*sin(2*ha*DEG2RAD);
                Dci=(12/PI)*(Fdi + Fdiss)*Dcz;
              }
              // Sum them up
              Gc=Bc+Dc;
              Gi=Bci+Dci;
            }

            if (beam_name) {
              if (is_not_null)
                beambuf[col] = Bci;
              else
                Rast_set_f_null_value(beambuf+col, 1);
            }

            if (diffuse_name) {
              if (is_not_null)
                diffusebuf[col] = Dci;
              else
                Rast_set_f_null_value(diffusebuf+col, 1);
            }

            if (total_name) {
              if (is_not_null)
                totalbuf[col] = Gi;
              else
                Rast_set_f_null_value(totalbuf+col, 1);
            }

            if (sunhour_name) {
              if (is_not_null) {
                sunhourbuf[col] = (pd.ssetr - pd.sretr) / 60.;
                if (sunhourbuf[col] > 24.)
                  sunhourbuf[col] = 24.;
                if (sunhourbuf[col] < 0.)
                  sunhourbuf[col] = 0.;
                // testtest sunhourbuf[col] = pd.ssha;
              }
              else
                Rast_set_f_null_value(sunhourbuf+col, 1);
            }
        }
        if (beam_name)
            Rast_put_f_row(beam_fd, beambuf);
        if (diffuse_name)
          Rast_put_f_row(diffuse_fd, diffusebuf);
        if (total_name)
          Rast_put_f_row(total_fd, totalbuf);
        if (sunhour_name)
            Rast_put_f_row(sunhour_fd, sunhourbuf);
    }
    G_percent(1, 1, 2);

    if (beam_name) {
        Rast_close(beam_fd);
        /* writing history file */
        Rast_short_history(beam_name, "raster", &hist);
        Rast_command_history(&hist);
        Rast_write_history(beam_name, &hist);
    }
    if (diffuse_name) {
      Rast_close(diffuse_fd);
      /* writing history file */
      Rast_short_history(diffuse_name, "raster", &hist);
      Rast_command_history(&hist);
      Rast_write_history(diffuse_name, &hist);
    }
    if (total_name) {
      Rast_close(total_fd);
      /* writing history file */
      Rast_short_history(total_name, "raster", &hist);
      Rast_command_history(&hist);
      Rast_write_history(total_name, &hist);
    }
    if (sunhour_name) {
        Rast_close(sunhour_fd);
        /* writing history file */
        Rast_short_history(sunhour_name, "raster", &hist);
        Rast_command_history(&hist);
        Rast_write_history(sunhour_name, &hist);
    }

    G_done_msg(" ");

    exit(EXIT_SUCCESS);
}

void set_solpos_time(struct posdata *pdat, int year, int month, int day,
		     int hour, int minute, int second, int timezone)
{
    pdat->year = year;
    pdat->month = month;
    pdat->day = day;
    pdat->daynum = day;
    pdat->hour = hour;
    pdat->minute = minute;
    pdat->second = second;
    pdat->timezone = timezone;

    pdat->time_updated = 1;
    pdat->longitude_updated = 1;
}

void set_solpos_longitude(struct posdata *pdat, double longitude)
{
    pdat->longitude = longitude;

    pdat->longitude_updated = 1;
}

int roundoff(double *x)
{
    /* watch out for the roundoff errors */
    if (fabs(*x) > 1.0) {
        if (*x > 0.0)
            *x = 1.0;
        else
            *x = -1.0;

        return 1;
    }

    return 0;
}
