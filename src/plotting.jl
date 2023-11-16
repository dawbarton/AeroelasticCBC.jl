push!(PGFPlotsX.CUSTOM_PREAMBLE, "\\usepackage{unicode-math}")
push!(PGFPlotsX.CUSTOM_PREAMBLE, "\\setmainfont{TeX Gyre Heros}")
push!(PGFPlotsX.CUSTOM_PREAMBLE, "\\setmathfont{TeX Gyre Pagella Math}")
push!(PGFPlotsX.CLASS_OPTIONS, "12pt")

function fig_openloop()
    @pgf Axis(
        {
            xlabel="Flow velocity (m/s)",
            ylabel="Heave amplitude (mm)",
            width=135mm,
            height=60mm,
            xmin=12,
            ymin=0,
            enlarge_x_limits=false,
        },
        PlotInc(

        )
    )
#     \begin{figure}
#     \centering
#     \begin{tikzpicture}
#       \begin{axis}[
#         xlabel={Flow velocity (m/s)},
#         ylabel={Heave amplitude (mm)},
#         width=135mm,
#         height=60mm,
#         xmin=12,
#         ymin=0,
#         enlarge x limits=false,
#         font={\sffamily},
#       ]
#         \addplot [mark=*, color=black, mark size=1pt, only marks]
#           table [x=velocity, y=heave, col sep=comma] {./Data/openloop_lco_day1.csv};
#         \addplot [mark=*, color=black, mark size=1pt, only marks]
#           table [x=velocity, y=heave, col sep=comma] {./Data/openloop_lco_day2.csv};
#         \addplot [mark=*, color=black, mark size=1pt, only marks]
#           table [x=velocity, y=heave, col sep=comma] {./Data/openloop_lco_day3.csv};
#         \addplot [mark=*, color=black, mark size=1pt, only marks]
#           table [x=velocity, y=heave, col sep=comma] {./Data/openloop_lco_day4.csv};
#         \addplot [mark=*, color=black, mark size=1pt, only marks]
#           table [x=velocity, y=heave, col sep=comma] {./Data/openloop_lco_day5.csv};
#         \addplot [mark=*, color=black, mark size=1pt, only marks]
#           table [x=velocity, y=heave, col sep=comma] {./Data/openloop_lco_day6.csv};
#         \addplot [mark=*, color=black, mark size=1pt, only marks]
#           table [x=velocity, y=heave, col sep=comma] {./Data/openloop_eq.csv};
#       \end{axis}
#     \end{tikzpicture}
#     \caption{The response of the flutter rig in open-loop mode (no control active) measured across multiple days. Bistable behaviour is evident for flow velocities between approximately 14\,m/s and 25\,m/s where a stable limit cycle oscillation in heave and pitch coexists with a stable equilibrium.}
#     \label{fig:openloop}
#   \end{figure}
#   @pgf

end
