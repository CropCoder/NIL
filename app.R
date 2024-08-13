library(shiny)
library(shinyWidgets)
library(bslib)
library(DT)
library(tidyverse)
library(openxlsx)
library(RIdeogram)
library(shinyalert)
library(colourpicker)
library(shinyjs)


source("scr_NIL_plot.R")

ui <- page_fixed(title = "NIL Tool 近等基因系差异位点分析与可视化工具",
                 theme = bs_theme(bootswatch = "flatly",
                                  # base_font = font_google("Inter"),
                                  primary = "#08979c",
                                  success = "#00474f",
                                  navbar_bg = "#203864"),
                 useShinyjs(),
                 # useShinyalert(),
                 tags$head(
                   tags$style(HTML("
      .light-green-background {
        background-color: #d9f7be; /* 淡绿色背景 */
        padding: 10px;           /* 内边距 */
        border-radius: 5px;      /* 边角圆滑 */
        border: 0px solid #cccccc; /* 边框 */
      }
    "))
                 ),
  br(),
  card(
      #h1("NIL 近等基因系差异位点分析与可视化工具",style="text-align: center"),
      tags$img(src = "https://jewin.oss-cn-hangzhou.aliyuncs.com/image-20240813172307118.png", width = "auto")
  ),
  br(),
  card(
      downloadButton("downloadData", "[ 工具测试专用 ]下载示例数据文件（测序获得的近等基因系变异位点基因型信息）"),
      fileInput("file_xlsx","请先点击下方按钮上传基因型文件 (xlsx)",width = "100%",buttonLabel = "打开文件",
                placeholder = "提示：请上传xlsx格式文件,按照位置排序的基因型文件，并只保留一个sheet，默认读取xlsx的sheet1",
                accept = ".xlsx"),
      p("提示：点击上方上传文件后，系统自动预览，请确认数据无误后开始分析(前三列依次是SNP名称、染色体、物理位置、参考基因型，右侧相邻两列为一对NIL基因型)"),
      hr(),
      accordion_panel("请点击这里预览刚刚上传的数据是否正确",icon = icon("table"),
        card(
          DTOutput("df_xlsx"),
          hr(),
          textOutput("info_SNP"),
          textOutput("info_sample"),
          textOutput("info_check")
        )
      ),
      layout_column_wrap(
          1/2,
          card(
              card_header("分析控制台"),
              radioGroupButtons(
                inputId = "job_type",
                label = NULL,
                choices = c("没有双亲数据","有双亲数据"),
                justified = TRUE,
                checkIcon = list(
                  yes = icon("ok", 
                             lib = "glyphicon"))
              ),
              conditionalPanel(
                condition = "input.job_type == '没有双亲数据'",
                div(
                  class = "light-green-background",
                  p("提示：仅有一对近等基因系情况下，采用该模式分析A和B样本之间的差异位点")
                ),
                layout_column_wrap(
                  1/2,
                  selectInput("opt_sam_A","请选择A样本名称",choices = NULL,multiple = F),
                  selectInput("opt_sam_B","请选择B样本名称",choices = NULL,multiple = F)
                )
              ),
              conditionalPanel(
                condition = "input.job_type == '有双亲数据'",
                div(
                  class = "light-green-background",
                  p("(未完成)提示：提供亲本的基因型后，系统会先分析子代A和B之间的差异位点，并根据亲本基因型分析子代基因型来源于哪个亲本，分别用不同颜色进行表示")
                ),
                layout_column_wrap(
                  1/2,
                  selectInput("opt_qin_1","请选择亲本1名称",choices = NULL,multiple = F),
                  selectInput("opt_qin_2","请选择亲本2名称",choices = NULL,multiple = F)
                ),
                layout_column_wrap(
                  1/2,
                  selectInput("opt_qintype_sam_A","请选择子代A样本名称",choices = NULL,multiple = F),
                  selectInput("opt_qintype_sam_B","请选择子代B样本名称",choices = NULL,multiple = F)
                ),
                layout_column_wrap(
                  1/2,
                  colourpicker::colourInput("opt_color_qin_1", "亲本1的代表颜色", "#389e0d"),
                  colourpicker::colourInput("opt_color_qin_2", "亲本2的代表颜色", "#ff7a45")
                )
              ),
              layout_column_wrap(
                1/2,
                colourpicker::colourInput("opt_color_diff", "选择差异位点展示颜色", "#747d8c"),
                colourpicker::colourInput("opt_color_zahe", "选择杂合位点展示颜色", "#f0f0f0")
              ),
              div(
                class = "light-green-background",
                p("扩展宽度越大，图中线条越粗。滑窗区间越小，保留的线越少")
              ),
              layout_column_wrap(
                  1/2,
                  numericInput("opt_line_width","单个位点的扩展显示宽度（kb）",value = 1), # 10Kb
                  numericInput("opt_window","滑窗区间大小（Mb）",value = 10) # 10MB
              ),
              materialSwitch(
                inputId = "need_PDF",
                label = "默认采用SVG可编辑格式，是否需要转换PDF（更加耗时）", inline = T,
                value = F,
                status = "success"
              ),
              actionButton("run","点击此处开始运行",icon = icon("check")),
              
              p("注意：请先运行分析结束后再下载结果文件"),
              downloadButton("down_xlsx","下载结果表格Excel"),
              downloadButton("down_svg","推荐下载结果图SVG格式（可编辑）"),
              conditionalPanel(
                condition = "input.need_PDF == true",
                downloadButton("down_pdf","下载结果图PDF格式（不可以编辑）")
              )
          ),
          card(
              includeMarkdown("intro.Rmd")
          )
      ),

      card(
        card_header("分析结果"),
        markdown(textOutput("out_stat_info"))
      ),
      # tags$img(src = "chromosome.png", width = "100%"),
      card(full_screen = T,
           # height = 500,
          # plotOutput("out_svg"),
          DTOutput("sata_info"),height = 800
      ),
      div(
        class = "light-green-background",
        p("开发者：赵记稳 （邮箱zhaojiwen@nwafu.edu.cn） ，欢迎反馈使用体验与错误信息")
      ),
      tags$img(src = "https://jewin.oss-cn-hangzhou.aliyuncs.com/image-20240706184430230.png", width = "auto")
  )
)

server <- function(input, output, session) {

  # 使用downloadHandler来管理文件下载
  output$downloadData <- downloadHandler(
    filename = function() {
      # 定义下载的文件名，可以根据需要动态修改
      paste("ExampleData.xlsx")
    },
    content = function(file) {
      # 在这里生成文件内容
      file.copy("example.xlsx",file)
    }
  )

  
  observeEvent(input$file_xlsx, {
      req(input$file_xlsx)
      df <- read.xlsx(input$file_xlsx$datapath, sheet = 1)
      
      
      colnames_sam <- colnames(df)[-c(1:4)]
      
      output$df_xlsx <- renderDT({
          df[1:7,]
      })
      
      output$info_SNP <- renderText(
          paste0("总变异位点数量（SNP）：",nrow(df),"个，其中含有NA的位点数量为",sum(is.na(df[,5])),"个。")
      )
      output$info_sample <- renderText(
          paste0("材料清单：(",length(colnames_sam),")  ",str_c(colnames_sam,sep="、",collapse = "、"))
      )
      
      
      output$info_check <- renderText(
        paste0("提示：您上传的 ",input$file_xlsx$filename," 文件符合格式，检查没有发现问题，可以进行分析")
      )
      
      updateSelectInput(session, "opt_sam_A", choices = colnames_sam)
      updateSelectInput(session, "opt_sam_B", choices = colnames_sam)
      
      shinyalert("上传成功","系统已成功解析EXCEL文件","success")
  })
    
    observeEvent(input$run,{
        req(input$file_xlsx)
        shinyalert("正在运行", "请勿刷新或关闭页面，等待运行结束会弹出提示信息......", type = "info")
        df <- read.xlsx(input$file_xlsx$datapath,sheet = 1)
        
        sample_A = input$opt_sam_A
        sample_B = input$opt_sam_B

        if (sample_A == sample_B){
          shinyalert("不能选择相同的两个样本", "请修改样本A和样本B，不能输入相同样本编号", type = "warning")

        }else{
          outdata <- get_NIL_plot(df,
                                sample_A = input$opt_sam_A,
                                sample_B = input$opt_sam_B,
                                color_diff = input$opt_color_diff,
                                color_zahe = input$opt_color_zahe,
                                line_width = input$opt_line_width,
                                opt_window = input$opt_window,
                                need_PDF = input$need_PDF
                                )
        
        output$sata_info <- renderDT(
                                     datatable(outdata$df_out, extensions = 'Buttons',
                                     options = list(
                                       dom = 'Bfrtip',
                                       buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                       pageLength = 13, # 设置默认的行数
                                       autoWidth = TRUE, # 自动调整列宽
                                       searching = TRUE, # 允许搜索
                                       filter = 'top' # 在每列的顶部添加过滤器
                                     ))
                            )
        
        output$out_stat_info <- renderText(outdata$sata_info)
        
        output$down_xlsx <- downloadHandler(
            filename = function(){
                paste0("OUT_Type_",sample_A,"_",sample_B,".xlsx")
            },
            content = function(file){
                file.copy(paste0("OUT_Type_",sample_A,"_",sample_B,".xlsx"),file)
            }
        )
        
        output$down_svg <- downloadHandler(
            filename = function(){
                paste0("Plot_SVG_",sample_A,"_",sample_B,".svg")
            },
            content = function(file){
                file.copy(paste0("Plot_SVG_",sample_A,"_",sample_B,".svg"),file)
            }
        )
        
        output$down_pdf <- downloadHandler(
            filename = function(){
                paste0("Plot_PDF_",sample_A,"_",sample_B,".pdf")
            },
            content = function(file){
                file.copy(paste0("Plot_PDF_",sample_A,"_",sample_B,".pdf"),file)
            }
        )

        shinyalert("完成", "运行已完成，请查看并下载结果", type = "success")
        }
        
        
    })
  
}

shinyApp(ui, server)