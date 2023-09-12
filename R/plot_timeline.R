plot_timeline <- function(){
    data_raw <- fread("data/data_cases_hosps_deaths_serology.csv")
    pop <- fread("data/population.csv")
    
    p1 <- ggplot(data_raw,aes(x = covid_multi_strain_date_as_date(day))) +
        geom_col(aes(y = cases,fill = "Cases"), alpha = 0.3) +
        geom_line(aes(y = 10*hosps, color = "Hospitalisations")) +
        geom_line(aes(y = 10*deaths, color = "Deaths")) +
        scale_y_continuous(name = "Cases",sec.axis = sec_axis(~./10, name = "Hospitalisations/Deaths")) +
        scale_fill_manual(name = "", values = c("Cases" = "blue")) +
        scale_color_manual(name = "", values = c("Hospitalisations" = "red","Deaths" = "black"), breaks = c("Hospitalisations", "Deaths")) + 
        labs(x = "Date") +
        theme_cowplot(font_size = 12) +
        theme(#axis.line.y = element_line(color = "blue"),
            #axis.text.y.left = element_text(color = "blue"),
            #axis.title.y.left = element_text(color = "blue"),
            #axis.line.y.right = element_line(color = "red"),
            #axis.text.y.right = element_text(color = "red"),
            #axis.title.y.right = element_text(color = "red")
            legend.position = c(0.01,0.97),
            legend.direction = "horizontal",
            legend.box = "horizontal"
        )
    
    p1_annttd <- p1 + annotate("text", x = as.Date("2020-07-25"), y = 650, label = "Borders reopened to\ninternational tourism +\nall quarantine\nmeasures lifted") + 
        annotate(geom = "segment", x = as.Date("2020-07-15"), y = 400, xend = as.Date("2020-07-15"), yend = 100, arrow = arrow(length = unit(2,"mm"))) +
        annotate("text", x = as.Date("2020-10-24"), y = 1050, label = "Curfew started on\nTahiti and Moorea") + 
        annotate(geom = "segment", x = as.Date("2020-10-24"), y = 900, xend = as.Date("2020-10-24"), yend = 600, arrow = arrow(length = unit(2,"mm"))) +
        annotate("text", x = as.Date("2021-01-18"), y = 500, label = "Vaccination campaign\nstarted") + 
        annotate(geom = "segment", x = as.Date("2021-01-18"), y = 400, xend = as.Date("2021-01-18"), yend = 100, arrow = arrow(length = unit(2,"mm"))) +
        annotate("text", x = as.Date("2021-05-01"), y = 550, label = "Borders\nreopened\nto tourism") + 
        annotate(geom = "segment", x = as.Date("2021-05-01"), y = 400, xend = as.Date("2021-05-01"), yend = 100, arrow = arrow(length = unit(2,"mm"))) +
        annotate("text", x = as.Date("2021-06-01"), y = 1050, label = "Increase in\ninternational flights") + 
        annotate(geom = "segment", x = as.Date("2021-06-01"), y = 900, xend = as.Date("2021-06-01"), yend = 100, arrow = arrow(length = unit(2,"mm"))) +
        annotate("text", x = as.Date("2021-06-20"), y = 1550, label = "Curfew + state of health \nemergency started +\nSunday confinement\non Tahiti and Moorea") + 
        annotate(geom = "segment", x = as.Date("2021-08-02"), y = 1300, xend = as.Date("2021-08-02"), yend = 800, arrow = arrow(length = unit(2,"mm"))) +
        annotate("text", x = as.Date("2021-09-20"), y = 1600, label = "Booster rollout\naccelerated") + 
        annotate(geom = "segment", x = as.Date("2021-08-23"), y = 1500, xend = as.Date("2021-08-23"), yend = 1300, arrow = arrow(length = unit(2,"mm"))) +
        annotate("text", x = as.Date("2021-11-15"), y = 550, label = "Curfew + state of health \nemergency lifted") + 
        annotate(geom = "segment", x = as.Date("2021-11-15"), y = 400, xend = as.Date("2021-11-15"), yend = 100, arrow = arrow(length = unit(2,"mm"))) + 
        annotate(geom = "text",x = as.Date("2020-10-21"), y = 2450, label = "Wild-type") + 
        annotate(geom = "segment",x = as.Date("2020-07-13"), y = 2350, xend = as.Date("2021-07-01"), yend = 2350, arrow = arrow(ends = "both", length = unit(2,"mm"))) + 
        annotate(geom = "text",x = as.Date("2021-08-15"), y = 2500, label = "Delta") +
        annotate(geom = "segment",x = as.Date("2021-06-11"), y = 2400, xend = as.Date("2022-02-06"), yend = 2400, arrow = arrow(ends = "both", length = unit(2,"mm"))) + 
        annotate(geom = "text",x = as.Date("2022-02-01"), y = 2550, label = "Omicron BA.1/BA.2") +
        annotate(geom = "segment",x = as.Date("2021-11-22"), y = 2450, xend = as.Date("2022-05-06"), yend = 2450, arrow = arrow(ends = "first", length = unit(2,"mm")))
    
    pop <- pop[,.(population = sum(total)),by = .(age_group)]
    p2 <- ggplot(pop,aes(x = population, y = age_group)) + 
        geom_col() + 
        labs(x = "Population",y = "Age group (years)") + 
        theme_cowplot(font_size = 12) +
        theme(axis.text.x = element_text(angle = 45,hjust = 1))
    
    p3 <- ggdraw() +
        draw_plot(p1_annttd) +
        draw_plot(p2, x = 0.07, y = .53, width = .15, height = .3)
    
    return(p3)
}

# data_raw_long <- melt(data_raw,id.vars = "day",measure.vars = patterns(cases = "cases_",hosps = "hosps_",deaths = "deaths_"))
# data_raw_long[,age_group := age_groups[variable]]
# outcomes_by_age <- data_raw_long[,lapply(.SD,function(x) sum(x,na.rm = T)),.SDcols = c("cases","hosps","deaths"),by = .(age_group)]
# 
# plot_grid(p11, p2, rel_widths = c(1,0.4), labels = "auto")
