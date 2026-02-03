# Load data
load_data_Tiller <- function() {
        
    df <- list()
    for (i in 1:12) {
        csv_path <- paste0("data/Tiller_edit/series_", i, ".csv")
        df_i_raw <- read_csv(csv_path, show_col_types = FALSE)
        #print(sort(names(df_i_raw)))
        df_i <- with(df_i_raw, {data.frame(
            t = t,
            F = F_co_A2D_vol,
            cgina = cginA,
            Fgina = FginA,
            Q = P_reb,
            yga = ygA,
            Na = N_A,
            cga = cgA,
            ca = clA,
            cd = clD
        )
        })
        df[[i]] <- df_i
        names(df)[i] <- paste0("series_", i)
    }
    return(df)

}