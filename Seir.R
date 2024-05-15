#' ODE solution
#'
#' This function solves an ODE system based on given parameter list.
#' @param user created parameters
#' @return solution of an ODE function system
#' @export
#' @name model

# 设置中文显示支持
sysfonts::font_add("kaiti", "STKAITI.TTF")
showtext::showtext_auto()

# 定义模型函数
model <- function(t, y, param) {
  N <- param["N"]
  # 从参数列表中获得传染病五种状态（S、E、I1、I2 以及 R）
  S <- y[1]
  E <- y[2]
  I1 <- y[3] # 第一次传染
  I2 <- y[4] # 第二次传染
  R <- y[5]

  # 从参数列表中提取模型参数
  beta <- param["beta"] # 传染率
  mu <- param["mu"] # 人口自然死亡率
  gamma <- param["gamma"] # 康复率
  lamda <- param["lamda"] # 感染潜伏期
  delta <- param["delta"] # 再次感染率

  # 基于模型参数和状态信息创建传染病传播数学模型
  dSt <- mu * (N - S) - beta * S * (I1 + I2) / N
  dEt <- beta * S * (I1 + I2) / N - (mu + lamda) * E
  dI1t <- lamda * E - (mu + delta) * I1
  dI2t <- delta * I1 - mu * I2 - gamma * I2
  dRt <- gamma * I2 - mu * R

  # 汇总模型求解结果
  outcome <- c(dSt, dEt, dI1t, dI2t, dRt)

  # 返回传染病五种状态求解结果列表
  return(list(outcome))
}

#设置仿真参数
times = seq(0,156,by=1/7)
param <- c(N = 1, mu = 0.01, lamda = 5, beta = 0.5, gamma = 0.1, delta = 0.01)
init = c(S=0.9999,E=0.00008,I1=0.00002,I2=0.00001,R =0)

result <- deSolve::ode(y = init,
                       times = times,
                       func = model,
                       parms = param)
result <- as.data.frame(result)

tail(round(result,3.6),10)

#绘制迁移图
#'@export

# 定义每个状态的颜色
status_colors <- c("blue", "orange", "red", "green", "purple")

#绘制迁移图
seirplot <- ggplot2::ggplot(data = result, ggplot2::aes(x = time)) +
  ggplot2::geom_line(ggplot2::aes(y = S, color = "S"), linewidth = 1.5) +
  ggplot2::geom_line(ggplot2::aes(y = E, color = "E"), linewidth = 1.5) +
  ggplot2::geom_line(ggplot2::aes(y = I1, color = "I1"), linewidth = 1.5) +
  ggplot2::geom_line(ggplot2::aes(y = I2, color = "I2"), linewidth = 1.5) +
  ggplot2::geom_line(ggplot2::aes(y = R, color = "R"), linewidth = 1.5) +
  ggplot2::scale_color_manual(values = status_colors) +  # 使用定义的颜色向量
  ggplot2::labs(x = "时间", y = "比率", color = "状态") +
  ggplot2::ggtitle("疾病模型仿真")

#运行图形对象绘制图像结果
seirplot

# 保存为矢量图
ggplot2::ggsave(seirplot, file = "seir_plot.pdf", width = 10, height = 6,dpi = 300)
ggplot2::ggsave(seirplot, file = "seir_plot.svg", width = 10, height = 6,dpi = 300)
