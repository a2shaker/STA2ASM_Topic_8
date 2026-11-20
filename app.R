# STA2ASM Week 8: Bayesian inference explorer
# Shinylive-friendly: only shiny and Base R. No package installation at runtime.
library(shiny)

PAL <- list(teal="#236A6E", blue="#174F8A", rose="#C65D7B",
            purple="#8064A2", grey="#666666", pale="#EEF5F5")
fmt <- function(x, digits=3) formatC(x, format="f", digits=digits)
# Render mathematical notation in dynamic and static Shiny content.
mj <- function(s) withMathJax(HTML(s))

# Numerically integrate a curve and obtain its central posterior quantiles.
trap_area <- function(x, y) sum(diff(x) * (head(y,-1) + tail(y,-1))/2)
curve_quantile <- function(x, density, probs) {
  cumulative <- c(0, cumsum(diff(x)*(head(density,-1)+tail(density,-1))/2))
  cumulative <- cumulative / tail(cumulative,1)
  approx(cumulative, x, xout=probs, ties="ordered")$y
}

# Two-sided Bayesian one-sample t-test based on the noncentral t likelihood.
# The effect-size prior location can be changed independently of the null value.
bayes_t <- function(x, family="Cauchy", scale=1/sqrt(2), prior_df=3, mu_null=0, prior_location=0) {
  n <- length(x)
  t_obs <- (mean(x)-mu_null)/(sd(x)/sqrt(n))
  df <- n-1
  d <- seq(min(-12,prior_location-12*scale),max(12,prior_location+12*scale),length.out=3001)
  prior <- switch(family,
                  "Cauchy"=dcauchy(d,prior_location,scale),
                  "Normal"=dnorm(d,prior_location,scale),
                  "Student's t"=dt((d-prior_location)/scale,df=prior_df)/scale)
  # dt() is a density in the observed test statistic, evaluated for each delta.
  likelihood <- dt(t_obs,df=df,ncp=d*sqrt(n))
  product <- prior*likelihood
  marginal_h1 <- trap_area(d,product)
  posterior <- product/marginal_h1
  marginal_h0 <- dt(t_obs,df=df)
  bf10 <- marginal_h1/marginal_h0
  q <- curve_quantile(d,posterior,c(.025,.5,.975))
  negative <- trap_area(d[d<=0],posterior[d<=0])
  prior0 <- switch(family,
                   "Cauchy"=dcauchy(0,prior_location,scale),
                   "Normal"=dnorm(0,prior_location,scale),
                   "Student's t"=dt((0-prior_location)/scale,df=prior_df)/scale)
  post0 <- prior0*marginal_h0/marginal_h1
  list(delta=d,prior=prior,likelihood=likelihood,product=product,
       posterior=posterior,t=t_obs,bf10=bf10,prior0=prior0,post0=post0,
       q=q,negative=negative,n=n)
}

# Normal-normal model, with known observation SD, for the prior-sensitivity tab.
update_normal <- function(mu0, prior_sd, sample_mean, n, sigma=15) {
  v0 <- prior_sd^2
  v <- 1/(1/v0+n/sigma^2)
  c(mean=v*(mu0/v0+n*sample_mean/sigma^2),sd=sqrt(v))
}

# Draw a filled density curve without introducing a plotting package.
draw_curve <- function(x,y,col,heading,ylab="Density",xlab=expression(delta),
                       ylim=NULL,zero=TRUE) {
  if(is.null(ylim)) ylim <- c(0,max(y)*1.30)
  plot(x,y,type="n",main=heading,xlab=xlab,ylab=ylab,ylim=ylim,las=1,bty="l")
  polygon(c(x,rev(x)),c(y,rep(0,length(y))),col=adjustcolor(col,alpha.f=.13),border=NA)
  lines(x,y,col=col,lwd=2.5)
  if(zero) {
    abline(v=0,col="#555555",lty=2,lwd=1.8)
    legend("topright",legend=expression(H[0]*":"~delta==0),
           col="#555555",lty=2,lwd=1.8,bty="n",cex=0.83,
           inset=c(0.07,0.02),x.intersp=0.65)
  }
}

ui <- fluidPage(
  tags$head(tags$style(HTML('
    body {font-family:"Helvetica Neue",Arial,sans-serif;color:#333;background:#fff;}
    h1,h2,h3,h4 {color:#236A6E;font-weight:500;}
    .navbar-default {background:white;border:0;border-bottom:1px solid #dce6e6;border-radius:0;box-shadow:none;}
    .navbar-default .navbar-nav>.active>a,.navbar-default .navbar-nav>.active>a:hover {background:#EEF5F5;color:#236A6E;border-bottom:3px solid #236A6E;}
    .well {background:#EEF5F5;border:1px solid #cfdfe0;border-left:4px solid #236A6E;}
    .btn-primary {background:#236A6E;border-color:#236A6E;}
    .result {background:#EEF5F5;border-left:4px solid #236A6E;padding:12px 16px;margin:12px 0;}
    .note {color:#58696a;font-size:13px;margin:10px 0;}
    table.asm-table {width:100%;border-collapse:collapse;margin:12px 0;}
    table.asm-table th,table.asm-table td {padding:9px 12px;border-bottom:1px solid #dce6e6;text-align:left;}
    table.asm-table th {color:#236A6E;background:#EEF5F5;}
    .shiny-plot-output {max-width:100%;}
    .bf-highlight {border:1px solid #dce6e6;border-radius:7px;background:#fafcfc;padding:12px 16px;margin:12px 0 18px;}
    .bf-heading {font-weight:600;font-size:15px;margin-bottom:9px;color:#333;}
    .bf-label {font-size:12px;color:#667578;margin-bottom:3px;}
    .bf-number {font-size:21px;line-height:1.45;}
    .bf-conclusion {font-size:13px;color:#58696a;margin-top:8px;}
  '))),
  titlePanel("STA2ASM - Week 8: Bayesian Inference"),
  withMathJax(),
  p(class="note","Explore Bayesian updating, compare inference methods and investigate the influence of prior information."),
  navbarPage(title=NULL,id="tab",collapsible=TRUE,
             tabPanel(mj("Bayesian \\(t\\)-test explorer"),
                      fluidRow(
                        column(3,wellPanel(
                          h4("Generate observations"),
                          numericInput("pop_mean",mj("Population mean \\(\\mu\\)"),-8,step=1),
                          numericInput("pop_sd",mj("Population SD \\(\\sigma\\)"),15,min=.1,step=1),
                          numericInput("sample_n",mj("Sample size \\(n\\)"),30,min=3,max=500,step=1),
                          actionButton("new_sample","Generate new sample",class="btn-primary"),
                          tags$hr(),h4("Hypothesis test"),
                          numericInput("mu_null",mj("Hypothesised mean \\(\\mu_0\\)"),0,step=1),
                          uiOutput("hypotheses"),
                          tags$hr(),h4(mj("Prior for standardised effect \\(\\delta\\)")),
                          selectInput("prior_family",mj("Prior distribution for \\(\\delta\\)"),c("Cauchy","Normal","Student's t")),
                          numericInput("prior_location",mj("Prior centre \\(\\delta_{\\mathrm{prior}}\\)"),0,step=.1),
                          sliderInput("prior_scale",mj("Prior SD (Normal) / scale (Cauchy or Student's \\(t\\))"),
                                      min=.15,max=2.5,value=1/sqrt(2),step=.01),
                          conditionalPanel("input.prior_family == \"Student's t\"",
                                           numericInput("prior_df",mj("Prior degrees of freedom \\(\\nu\\)"),3,min=1,max=30,step=1)),
                          actionButton("restore_prior","Restore default Cauchy prior"),
                          tags$hr(),
                          p(class="note","Bootstrap inference uses B = 1,000 resamples for each generated sample.")
                        )),
                        column(9,
                               h3("Observed sample"),uiOutput("sample_summary"),
                               uiOutput("bf_highlight"),
                               h3("Visualising the Bayesian components"),
                               fluidRow(column(6,plotOutput("prior_plot",height="285px")),
                                        column(6,plotOutput("likelihood_plot",height="285px"))),
                               fluidRow(column(6,plotOutput("product_plot",height="285px")),
                                        column(6,plotOutput("posterior_plot",height="285px"))),
                               h3("Prior and posterior: Savage-Dickey"),
                               uiOutput("bayes_summary"),
                               fluidRow(
                                 column(4,checkboxInput("show_cri","Shade central 95% credible interval",TRUE)),
                                 column(4,checkboxInput("show_direction","Shade directional posterior probability",FALSE)),
                                 column(4,radioButtons("direction",mj("Direction for \\(\\delta\\)"),
                                                       choices=c("Below zero"="negative","Above zero"="positive"),
                                                       selected="negative",inline=TRUE))
                               ),
                               plotOutput("combined_plot",height="430px"),
                               h3("Comparison of inference methods"),
                               uiOutput("comparison")
                        )
                      )
             ),
             tabPanel("Prior sensitivity and sample size",
                      fluidRow(
                        column(3,wellPanel(
                          h4("Prior sensitivity"),
                          selectInput("scenario","Prior presets",choices=c(
                            "Broad prior"="broad","Moderately informative"="moderate",
                            "Informative, conflicting"="conflict","Strongly informative, conflicting"="strong",
                            "Custom"="custom")),
                          p(class="note", "Choose a preset to set the prior mean and SD, or adjust them to create a custom prior."),
                          numericInput("prior_mu",mj("Prior mean \\(\\mu_0\\)"),0,step=1),
                          sliderInput("prior_sd",mj("Prior SD \\(\\sigma_0\\)"),min=1,max=30,value=20,step=.5),
                          numericInput("observed_mu",mj("Observed sample mean \\(\\bar{x}\\)"),-9,step=1),
                          numericInput("known_sigma",mj("Assumed known population SD \\(\\sigma\\)"),15,min=.1,step=1),
                          sliderInput("n_small",mj("Small sample size \\(n_1\\)"),min=3,max=40,value=5,step=1),
                          sliderInput("n_large",mj("Large sample size \\(n_2\\)"),min=40,max=300,value=100,step=5),
                          mj("<p class='note'>The observed sample mean \\(\\bar{x}\\) is held constant across sample sizes, isolating the influence of the prior and sample size.</p>")
                        )),
                        column(9,
                               h3("How prior information changes the posterior"),
                               mj("<p>This tab uses a normal-normal model with known population standard deviation \\(\\sigma\\), as in the readings. It is different from the noncentral \\(t\\) model in Tab 1.</p>"),
                               plotOutput("sensitivity_plot",height="475px"),
                               uiOutput("sensitivity_summary"),
                               mj("<p class='note'>The posterior combines the prior and the observed data. A strongly informative prior can dominate a small sample, even when the prior conflicts with the observed sample mean \\(\\bar{x}\\).</p>")
                        )
                      )
             )
  )
)

server <- function(input,output,session) {
  # The initial simulated dataset retains the Week 8 numerical defaults.
  sample_data <- eventReactive(list(input$new_sample,input$pop_mean,input$pop_sd,input$sample_n),{
    if(input$new_sample==0 && input$pop_mean==-8 && input$pop_sd==15 && input$sample_n==30)
      set.seed(512)
    rnorm(as.integer(input$sample_n),mean=input$pop_mean,sd=input$pop_sd)
  },ignoreNULL=FALSE)
  
  observeEvent(input$restore_prior,{
    updateSelectInput(session,"prior_family",selected="Cauchy")
    updateSliderInput(session,"prior_scale",value=1/sqrt(2))
    updateNumericInput(session,"prior_df",value=3)
    updateNumericInput(session,"prior_location",value=0)
  })
  
  analysis <- reactive({
    req(sample_data(),input$prior_family,input$prior_scale,input$prior_location)
    bayes_t(sample_data(),input$prior_family,input$prior_scale,
            if(is.null(input$prior_df)) 3 else input$prior_df,
            mu_null=input$mu_null,prior_location=input$prior_location)
  })
  
  boot_result <- eventReactive(list(sample_data(),input$mu_null),{
    x <- sample_data(); n <- length(x); B <- 1000
    # Ordinary percentile bootstrap CI for the population mean.
    means <- replicate(B,mean(sample(x,n,replace=TRUE)))
    ci <- unname(quantile(means,c(.025,.975)))
    # Null-centred bootstrap test, studentised on each resample.
    x0 <- x-mean(x)+input$mu_null
    t0 <- replicate(B,{
      xb <- sample(x0,n,replace=TRUE)
      (mean(xb)-input$mu_null)/(sd(xb)/sqrt(n))
    })
    tobs <- (mean(x)-input$mu_null)/(sd(x)/sqrt(n))
    p <- (1+sum(abs(t0)>=abs(tobs)))/(B+1)
    list(ci=ci,p=p)
  },ignoreNULL=FALSE)
  
  output$hypotheses <- renderUI({
    req(input$mu_null)
    mj(sprintf(paste0(
      "<div class='note'><strong>Two-sided hypotheses</strong><br>",
      "\\(H_0:\\mu=%s\\) &nbsp; vs &nbsp; \\(H_1:\\mu\\ne %s\\)<br>",
      "Equivalently, for \\(\\delta=(\\mu-\\mu_0)/\\sigma\\):<br>",
      "\\(H_0:\\delta=0\\) &nbsp; vs &nbsp; \\(H_1:\\delta\\ne0\\)</div>"),
      fmt(input$mu_null),fmt(input$mu_null)))
  })
  
  output$sample_summary <- renderUI({
    x <- sample_data()
    mj(sprintf(
      "<div class='result'>\\(n=%d\\) &nbsp; | &nbsp; Sample mean: \\(\\bar{x}=%s\\) &nbsp; | &nbsp; Sample SD: \\(s=%s\\) &nbsp; | &nbsp; \\(t_{%d}=%s\\)</div>",
      length(x),fmt(mean(x)),fmt(sd(x)),length(x)-1,fmt(analysis()$t)))
  })
  
  output$prior_plot <- renderPlot({
    a<-analysis(); draw_curve(a$delta,a$prior,PAL$blue,"1. Prior")
  })
  output$likelihood_plot <- renderPlot({
    a<-analysis(); draw_curve(a$delta,a$likelihood,PAL$rose,"2. Likelihood","Likelihood")
  })
  output$product_plot <- renderPlot({
    a<-analysis(); draw_curve(a$delta,a$product,PAL$grey,"3. Prior x likelihood","Unscaled density")
  })
  # The four components use a common scale to show how updating works.
  # Posterior annotations and shaded probabilities are reserved for the larger plot.
  output$posterior_plot <- renderPlot({
    a <- analysis()
    draw_curve(a$delta,a$posterior,PAL$teal,"4. Posterior")
  })
  
  output$combined_plot <- renderPlot({
    a <- analysis()
    d <- a$delta
    
    # Adaptive limits: include the central 84% of the prior, almost all of
    # the posterior, and delta = 0 for the Savage-Dickey markers.
    # This avoids displaying the uninformative extreme tails of a Cauchy prior.
    prior_bounds <- switch(input$prior_family,
                           "Cauchy"=qcauchy(c(.08,.92),location=input$prior_location,scale=input$prior_scale),
                           "Normal"=qnorm(c(.08,.92),mean=input$prior_location,sd=input$prior_scale),
                           "Student's t"=input$prior_location + input$prior_scale *
                             qt(c(.08,.92),df=if(is.null(input$prior_df)) 3 else input$prior_df))
    post_bounds <- a$q[c(1,3)]
    # Extend the posterior interval to cover almost all its mass.
    post_bounds <- curve_quantile(d,a$posterior,c(.001,.999))
    xr <- range(c(0,prior_bounds,post_bounds))
    pad <- max(diff(xr)*.10,0.12)
    xlim <- xr + c(-pad,pad)
    
    visible <- d >= xlim[1] & d <= xlim[2]
    ylim <- c(0,max(c(a$prior[visible],a$posterior[visible],
                      a$prior0,a$post0))*1.22)
    
    plot(d,a$prior,type="n",xlim=xlim,ylim=ylim,
         xlab=expression("Standardised effect "*delta),ylab="Probability density",
         main="Prior and posterior under H1",las=1,bty="l")
    
    # Shade the central 95% credible interval (equal-tailed).
    if(isTRUE(input$show_cri)) {
      xs <- sort(unique(c(a$q[1],d[d>a$q[1] & d<a$q[3]],a$q[3])))
      ys <- approx(d,a$posterior,xout=xs)$y
      polygon(c(xs[1],xs,tail(xs,1)),c(0,ys,0),
              col=adjustcolor(PAL$teal,alpha.f=.23),border=NA)
      abline(v=a$q[c(1,3)],col=PAL$rose,lty=3,lwd=1.7)
    }
    
    # Shade P(delta < 0 | data, H1) or P(delta > 0 | data, H1).
    if(isTRUE(input$show_direction)) {
      bounds <- if(input$direction=="positive") c(0,max(d)) else c(min(d),0)
      xs <- sort(unique(c(bounds[1],d[d>bounds[1] & d<bounds[2]],bounds[2])))
      ys <- approx(d,a$posterior,xout=xs)$y
      polygon(c(xs[1],xs,tail(xs,1)),c(0,ys,0),
              col=adjustcolor(PAL$purple,alpha.f=.27),border=NA)
    }
    
    lines(d,a$prior,col=PAL$blue,lwd=2.7,lty=2)
    lines(d,a$posterior,col=PAL$teal,lwd=2.7)
    abline(v=0,col="#555555",lty=2,lwd=1.8)
    abline(v=a$q[2],col=PAL$teal,lty=3,lwd=1.5)
    points(0,a$prior0,pch=21,bg="grey70",col="grey35",cex=1.6)
    points(0,a$post0,pch=21,bg="grey70",col="grey35",cex=1.6)
    
    legend("topright",
           legend=c("Prior","Posterior",expression(H[0]*":"~delta==0),
                    "Density at delta = 0","Posterior median",
                    if(isTRUE(input$show_cri)) "95% CrI" else NULL,
                    if(isTRUE(input$show_direction)) "Directional probability" else NULL),
           col=c(PAL$blue,PAL$teal,"#555555","grey35",PAL$teal,
                 if(isTRUE(input$show_cri)) PAL$rose else NULL,
                 if(isTRUE(input$show_direction)) PAL$purple else NULL),
           lty=c(2,1,2,NA,3,
                 if(isTRUE(input$show_cri)) 3 else NULL,
                 if(isTRUE(input$show_direction)) 1 else NULL),
           lwd=c(2.7,2.7,1.8,NA,1.5,
                 if(isTRUE(input$show_cri)) 1.7 else NULL,
                 if(isTRUE(input$show_direction)) 4 else NULL),
           pch=c(NA,NA,NA,21,NA,
                 if(isTRUE(input$show_cri)) NA else NULL,
                 if(isTRUE(input$show_direction)) NA else NULL),
           pt.bg=c(NA,NA,NA,"grey70",NA,
                   if(isTRUE(input$show_cri)) NA else NULL,
                   if(isTRUE(input$show_direction)) NA else NULL),
           bty="n",cex=.85,inset=c(0.025,0.02),x.intersp=0.7)
  })
  output$bf_highlight <- renderUI({
    a <- analysis()
    bf <- a$bf10
    if (!is.finite(bf) || bf <= 0) return(NULL)
    preferred <- if (bf >= 1) "H_1" else "H_0"
    main_bf <- if (bf >= 1) "\\text{BF}_{10}" else "\\text{BF}_{01}"
    other_bf <- if (bf >= 1) "\\text{BF}_{01}" else "\\text{BF}_{10}"
    main_value <- max(bf, 1/bf)
    other_value <- min(bf, 1/bf)
    strength <- if (main_value > 100) "Decisive" else if (main_value >= 30) "Very strong" else if (main_value >= 10) "Strong" else if (main_value >= 3) "Moderate" else "Anecdotal"
    message <- if (abs(bf-1) < 1e-8) "Equal evidence for both hypotheses." else
      sprintf("%s evidence in favour of \\(%s\\).", strength, preferred)
    div(class="bf-highlight",
        tags$div(class="bf-heading", "Bayesian hypothesis test"),
        fluidRow(
          column(6,
                 tags$div(class="bf-label", if (bf >= 1) "Evidence favouring H1" else "Evidence favouring H0"),
                 mj(sprintf("<div class='bf-number'>\\(%s=%s\\)</div>", main_bf, fmt(main_value)))
          ),
          column(6,
                 tags$div(class="bf-label", "Reciprocal"),
                 mj(sprintf("<div class='bf-number'>\\(%s=%s\\)</div>", other_bf, fmt(other_value)))
          )
        ),
        mj(sprintf("<div class='bf-conclusion'>%s</div>", message))
    )
  })
  output$bayes_summary <- renderUI({
    a<-analysis(); bf<-a$bf10; probh1<-bf/(1+bf)
    direction_prob <- if(input$direction=="positive") 1-a$negative else a$negative
    direction_sign <- if(input$direction=="positive") ">" else "<"
    div(class="result",
        tags$strong("Bayesian hypothesis testing"),tags$br(),
        withMathJax(HTML(sprintf(
          "\\(\\text{BF}_{10}=%s\\), \\(\\text{BF}_{01}=%s\\). With prior odds 1:1, \\(P(H_1\\mid\\mathrm{data})=%s\\) and \\(P(H_0\\mid\\mathrm{data})=%s\\).",
          fmt(bf),fmt(1/bf),fmt(probh1),fmt(1-probh1)))),
        tags$br(),tags$br(),tags$strong("Savage-Dickey density ratio"),tags$br(),
        mj(sprintf("Prior density at \\(\\delta=0\\) = %s; posterior density at \\(\\delta=0\\) = %s; their ratio = %s (\\(\\text{BF}_{10}\\)).",
                   fmt(a$prior0,5),fmt(a$post0,5),fmt(a$prior0/a$post0,3))),
        tags$br(),tags$br(),mj("<strong>Posterior estimation, conditional on \\(H_1\\)</strong>"),tags$br(),
        withMathJax(HTML(sprintf(
          "Median of \\(\\delta\\): %s. Central 95%% credible interval: (%s, %s). \\(P(\\delta%s0\\mid\\mathrm{data},H_1)=%s\\).",
          fmt(a$q[2]),fmt(a$q[1]),fmt(a$q[3]),direction_sign,fmt(direction_prob,4))))
    )
  })
  output$comparison <- renderUI({
    x <- sample_data(); tr <- t.test(x,mu=input$mu_null); br <- boot_result(); a <- analysis()
    prob_h1 <- a$bf10/(1+a$bf10)
    # Display the Bayes factor favouring the better-supported model first.
    preferred_bf <- if (a$bf10 >= 1) {
      sprintf("\\text{BF}_{10}=%s", fmt(a$bf10))
    } else {
      sprintf("\\text{BF}_{01}=%s", fmt(1/a$bf10))
    }
    reciprocal_bf <- if (a$bf10 >= 1) {
      sprintf("\\text{BF}_{01}=%s", fmt(1/a$bf10))
    } else {
      sprintf("\\text{BF}_{10}=%s", fmt(a$bf10))
    }
    strength_value <- max(a$bf10, 1/a$bf10)
    strength <- if (strength_value > 100) "Decisive" else if (strength_value >= 30) "Very strong" else if (strength_value >= 10) "Strong" else if (strength_value >= 3) "Moderate" else "Anecdotal"
    preferred_h <- if (a$bf10 >= 1) "H_1" else "H_0"
    evidence_label <- if (abs(a$bf10 - 1) < 1e-8) "Equal evidence for both hypotheses" else sprintf("%s evidence in favour of \\(%s\\)", strength, preferred_h)
    bf_text <- sprintf("\\(%s;\\ %s\\). %s", preferred_bf, reciprocal_bf, evidence_label)
    directional <- if(input$direction=="positive") 1-a$negative else a$negative
    dir_sign <- if(input$direction=="positive") ">" else "<"
    # Under the normal sampling model with the reference non-informative
    # location-scale prior, the 95% CrI for mu equals the t confidence interval.
    mu_cri <- tr$conf.int
    rows <- list(
      c("Frequentist \\(t\\)-test",
        paste0("\\(p=",format.pval(tr$p.value,digits=3),"\\)"),
        sprintf("(%.2f, %.2f) for \\(\\mu\\)",tr$conf.int[1],tr$conf.int[2])),
      c("Bootstrap",
        paste0("\\(p=",format.pval(br$p,digits=3),"\\)"),
        sprintf("(%.2f, %.2f) for \\(\\mu\\)",br$ci[1],br$ci[2])),
      c("Bayesian \\(t\\)-test: effect size",
        sprintf("%s; \\(P(H_1\\mid\\mathrm{data})=%s\\); \\(P(\\delta%s0\\mid\\mathrm{data},H_1)=%s\\)",
                bf_text,fmt(prob_h1,3),dir_sign,fmt(directional,4)),
        sprintf("(%.3f, %.3f) for \\(\\delta\\)",a$q[1],a$q[3])),
      c("Bayesian: population mean",
        "Non-informative prior; same numerical interval as frequentist \\(t\\)-test CI",
        sprintf("(%.2f, %.2f) for \\(\\mu\\)",mu_cri[1],mu_cri[2]))
    )
    mj(paste0(
      "<table class='asm-table'><thead><tr><th>Method</th><th>Hypothesis test / posterior probability</th><th>95% interval</th></tr></thead><tbody>",
      paste(vapply(rows,function(r) paste0("<tr>",paste0("<td>",r,"</td>",collapse=""),"</tr>"),character(1)),collapse=""),
      "</tbody></table>"))
  })
  
  # Presets populate the controls; manually editing either control selects Custom.
  prior_presets <- list(broad=c(0,20), moderate=c(0,8),
                        conflict=c(8,5), strong=c(15,3))
  observeEvent(input$scenario, {
    if (identical(input$scenario, "custom")) return()
    values <- prior_presets[[input$scenario]]
    req(values)
    freezeReactiveValue(input, "prior_mu")
    freezeReactiveValue(input, "prior_sd")
    updateNumericInput(session, "prior_mu", value=values[1])
    updateSliderInput(session, "prior_sd", value=values[2])
  }, ignoreInit=TRUE)
  
  observeEvent(list(input$prior_mu, input$prior_sd), {
    selected <- input$scenario
    if (is.null(selected) || identical(selected, "custom")) return()
    expected <- prior_presets[[selected]]
    if (!is.null(expected) &&
        (!isTRUE(all.equal(as.numeric(input$prior_mu), expected[1])) ||
         !isTRUE(all.equal(as.numeric(input$prior_sd), expected[2])))) {
      updateSelectInput(session, "scenario", selected="custom")
    }
  }, ignoreInit=TRUE)
  
  sensitivity <- reactive({
    req(input$prior_mu,input$prior_sd,input$observed_mu,input$known_sigma,
        input$n_small,input$n_large)
    small<-update_normal(input$prior_mu,input$prior_sd,input$observed_mu,
                         input$n_small,input$known_sigma)
    large<-update_normal(input$prior_mu,input$prior_sd,input$observed_mu,
                         input$n_large,input$known_sigma)
    list(small=small,large=large)
  })
  output$sensitivity_plot <- renderPlot({
    s<-sensitivity()
    centres<-c(input$prior_mu,input$observed_mu,s$small[1],s$large[1])
    spreads<-c(input$prior_sd,input$known_sigma/sqrt(input$n_small),
               s$small[2],s$large[2])
    lim<-range(centres+c(-1,1)*max(spreads)*3)
    d<-seq(lim[1],lim[2],length.out=1400)
    prior<-dnorm(d,input$prior_mu,input$prior_sd)
    psmall<-dnorm(d,s$small[1],s$small[2])
    plarge<-dnorm(d,s$large[1],s$large[2])
    plot(d,prior,type="l",col=PAL$blue,lwd=2.7,lty=2,las=1,bty="l",
         ylim=c(0,max(prior,psmall,plarge)*1.35),
         xlab="Population mean",ylab="Density",
         main="Prior and posterior distributions")
    lines(d,psmall,col=PAL$rose,lwd=2.7)
    lines(d,plarge,col=PAL$teal,lwd=2.7)
    abline(v=input$observed_mu,col="grey30",lty=2,lwd=1.8)
    legend("topright",legend=c("Prior",paste0("Posterior: n = ",input$n_small),
                               paste0("Posterior: n = ",input$n_large),
                               paste0("Observed mean = ",fmt(input$observed_mu,1))),
           col=c(PAL$blue,PAL$rose,PAL$teal,"grey30"),
           lty=c(2,1,1,2),lwd=c(2.7,2.7,2.7,1.8),bty="n",cex=.9)
  })
  output$sensitivity_summary <- renderUI({
    s <- sensitivity()
    ci <- function(z) sprintf("(%.2f, %.2f)",z[1]-qnorm(.975)*z[2],
                              z[1]+qnorm(.975)*z[2])
    mj(sprintf(paste0(
      "<div class='result'><strong>Posterior summaries</strong>",
      "<p class='note'>The mean below is the <em>posterior mean</em> \\(\\text{E}(\\mu\\mid\\mathrm{data})\\),",
      " our estimate of the population mean \\(\\mu\\), not the observed sample mean \\(\\bar{x}\\).</p>",
      "Small sample (\\(n=%d\\)): posterior mean \\(\\text{E}(\\mu\\mid\\mathrm{data})=%s\\); 95%% CrI = %s<br>",
      "Large sample (\\(n=%d\\)): posterior mean \\(\\text{E}(\\mu\\mid\\mathrm{data})=%s\\); 95%% CrI = %s",
      "</div>"),
      input$n_small,fmt(s$small[1]),ci(s$small),
      input$n_large,fmt(s$large[1]),ci(s$large)))
  })
  
}

shinyApp(ui,server)
