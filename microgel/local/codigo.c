    if (ludwig->hydro) {

      /* Zero velocity field here, as velocity at collision is used
       * at next time step for FD above. Strictly, we only need to
       * do this if velocity output is required in presence of
       * colloids to present non-zero u inside particles. */

      hydro_u_zero(ludwig->hydro, uzero);

      /* Collision stage */
      lb_collide(ludwig->lb, ludwig->hydro, ludwig->map, ludwig->noise,
     ludwig->fe, ludwig->visc);

      /* Colloid bounce-back applied between collision and
       * propagation steps. */
      wall_set_wall_distributions(ludwig->wall);
      subgrid_update(ludwig->collinfo, ludwig->hydro, ludwig->lb->param->noise);
      bounce_back_on_links(ludwig->bbl, ludwig->lb, ludwig->wall,
         ludwig->collinfo);
      wall_bbl(ludwig->wall);

    }

    /* There must be no halo updates between bounce back
     * and propagation, as the halo regions are active */

    if (ludwig->hydro) {
      lb_propagation(ludwig->lb);
    }

// Estadisticas