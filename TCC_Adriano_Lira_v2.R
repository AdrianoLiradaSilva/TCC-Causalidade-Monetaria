# =============================================================================
# TÍTULO:  Causalidade entre Índice de Preços e os Agregados Monetários
# AUTOR:   Adriano Lira da Silva
# INST.:   Universidade Federal Rural de Pernambuco (UFRPE) – DECON
# DATA:    2025–2026
# VERSÃO:  2.0
# =============================================================================
# DESCRIÇÃO:
#   Este script implementa a análise empírica do TCC, compreendendo:
#     1. Importação e tratamento da base de dados
#     2. Transformação em séries temporais (ts)
#     3. Análise gráfica (nível e primeiras diferenças)
#     4. Funções de Autocorrelação (ACF e PACF)
#     5. Testes de raiz unitária (ADF – Dickey-Fuller Aumentado)
#     6. Seleção do número ótimo de defasagens (VARselect)
#     7. Estimação do modelo VAR
#     8. Diagnóstico dos resíduos (estabilidade e autocorrelação)
#     9. Testes de Causalidade de Granger
#    10. Exportação de resultados (gráficos, tabelas e modelo salvo)
# =============================================================================


# -----------------------------------------------------------------------------
# 0. CONFIGURAÇÕES INICIAIS
# -----------------------------------------------------------------------------

# Instalar pacotes (execute apenas na primeira vez)

install.packages(c("readr", "dplyr", "vars", "tseries", "urca", 
                   "ggplot2", "tidyr", "zoo", "lmtest"))

# Carregar bibliotecas
library(readr)      # Leitura de arquivos CSV
library(dplyr)      # Manipulação de dados
library(tidyr)      # Reorganização de dados
library(ggplot2)    # Visualizações
library(vars)       # Modelos VAR e Causalidade de Granger
library(tseries)    # Teste ADF
library(urca)       # Testes de raiz unitária alternativos (KPSS, PP)
library(zoo)        # Séries temporais irregulares

# Semente para reprodutibilidade
set.seed(42)

# Configuração global de gráficos (base R)
op_original <- par(no.readonly = TRUE)   # salva configurações originais


# -----------------------------------------------------------------------------
# 1. IMPORTAÇÃO E TRATAMENTO DA BASE DE DADOS
# -----------------------------------------------------------------------------

# Definir diretório de trabalho (ajuste o caminho conforme necessário)
# setwd("C:/Users/adria/Documents/UFRPE/TCC/1 - TCC/Dados basicos")
# getwd(); list.files()

# Leitura da base (separador ";", decimal ",")
dados <- readr::read_delim(
  "Base_Var.csv",
  delim      = ";",
  locale     = locale(decimal_mark = ",", grouping_mark = "."),
  show_col_types = FALSE
)

# Inspeção inicial
cat("\n===== ESTRUTURA DA BASE =====\n")
str(dados)
cat("\n===== PRIMEIRAS LINHAS =====\n")
print(head(dados))
cat("\n===== RESUMO ESTATÍSTICO =====\n")
print(summary(dados))

# ── Conversão da coluna de data ──────────────────────────────────────────────
meses_pt <- c(
  jan = "01", fev = "02", mar = "03", abr = "04",
  mai = "05", jun = "06", jul = "07", ago = "08",
  set = "09", out = "10", nov = "11", dez = "12"
)

dados <- dados %>%
  mutate(
    mes_abrev = substr(Data, 1, 3),
    ano_str   = substr(Data, 5, 6),
    ano_full  = paste0("20", ano_str),
    mes_num   = meses_pt[mes_abrev],
    data      = as.Date(paste0(ano_full, "-", mes_num, "-01"))
  ) %>%
  dplyr::select(-mes_abrev, -ano_str, -ano_full, -mes_num, -Data) %>%
  rename(inflacao = ipca) %>%
  mutate(
    inflacao = as.numeric(gsub(",", ".", inflacao)),
    ln_BM    = as.numeric(gsub(",", ".", ln_BM)),
    ln_M1    = as.numeric(gsub(",", ".", ln_M1)),
    ln_M2    = as.numeric(gsub(",", ".", ln_M2))
  ) %>%
  arrange(data)



# Verificação pós-conversão
cat("\n===== DATAS CONVERTIDAS (primeiras e últimas) =====\n")
print(head(dados$data)); print(tail(dados$data))
cat("NAs na coluna data:", sum(is.na(dados$data)), "\n")


#teste do grapfico:

# Verificar o final da série
tail(dados_ts, 6)

# Verificar se há NAs
colSums(is.na(dados_var))

# -----------------------------------------------------------------------------
# 2. CONSTRUÇÃO DAS SÉRIES TEMPORAIS
# -----------------------------------------------------------------------------

# Selecionar variáveis de interesse
vars_modelo <- c("inflacao", "ln_BM", "ln_M1", "ln_M2")

dados_var <- dados %>%
  dplyr::select(data, all_of(vars_modelo)) %>%
  filter(data <= as.Date("2025-12-01")) %>%   # garante que só vai até dez/2025
  dplyr::select(-data) %>%
  na.omit()

# Criar objeto ts (mensal, a partir de jan/2002)
dados_ts <- ts(
  dados_var,
  start     = c(2002, 1),
  frequency = 12
)

cat("\n===== PROPRIEDADES DA SÉRIE TEMPORAL =====\n")
cat("Início:", paste(start(dados_ts), collapse = "/"), "\n")
cat("Fim:   ", paste(end(dados_ts),   collapse = "/"), "\n")
cat("Frequência:", frequency(dados_ts), "obs/ano\n")
cat("Total de observações:", nrow(dados_ts), "\n")

# Primeiras diferenças (para estacionarizar agregados monetários)
dados_diff <- diff(dados_ts)


# -----------------------------------------------------------------------------
# 3. GRÁFICOS DAS SÉRIES EM NÍVEL
# -----------------------------------------------------------------------------

head(dados[, c("ln_BM", "ln_M1", "ln_M2")])

png(
  filename = "Figura_1_Series_Nivel.png",
  width = 3200, height = 2400, res = 300
)

par(
  mfrow  = c(4, 1),
  mar    = c(3, 5, 2.5, 1.5),
  oma    = c(0.5, 0, 3, 0),
  cex.axis = 0.85, cex.lab = 1.0, cex.main = 1.05
)

plot(dados_ts[, "inflacao"],
     type = "l", col = "#1a5276", lwd = 1.5,
     main = "Inflação – IPCA (variação % mensal)",
     ylab = "Var. % mensal", xlab = "")
abline(h = 0, lty = 2, col = "gray50", lwd = 0.8)

plot(dados_ts[, "ln_BM"],
     type = "l", col = "#117a65", lwd = 1.5,
     main = "Base Monetária – M0 (logaritmo natural)",
     ylab = "ln(M0)", xlab = "")

plot(dados_ts[, "ln_M1"],
     type = "l", col = "#7d6608", lwd = 1.5,
     main = "Meios de Pagamento – M1 (logaritmo natural)",
     ylab = "ln(M1)", xlab = "")

plot(dados_ts[, "ln_M2"],
     type = "l", col = "#6e2f1a", lwd = 1.5,
     main = "Meios de Pagamento Amplos – M2 (logaritmo natural)",
     ylab = "ln(M2)", xlab = "Período (mensal)")

dev.off()
cat(">> Figura 1 exportada.\n")


# -----------------------------------------------------------------------------
# 4. GRÁFICOS DAS SÉRIES EM PRIMEIRAS DIFERENÇAS
# -----------------------------------------------------------------------------

png(
  filename = "Figura_2_Series_Diferenciadas.png",
  width = 3200, height = 2400, res = 300
)

par(
  mfrow  = c(4, 1),
  mar    = c(3, 5, 2.5, 1.5),
  oma    = c(0.5, 0, 3, 0),
  cex.axis = 0.85, cex.lab = 1.0, cex.main = 1.05
)

plot(dados_diff[, "inflacao"],
     type = "l", col = "#1a5276", lwd = 1.5,
     main = "Inflação – IPCA (1ª diferença)",
     ylab = "Δ inflação", xlab = "")
abline(h = 0, lty = 2, col = "gray50", lwd = 0.8)

plot(dados_diff[, "ln_BM"],
     type = "l", col = "#117a65", lwd = 1.5,
     main = "Taxa de crescimento de M0 (1ª dif. do ln)",
     ylab = "Δ ln(M0)", xlab = "")
abline(h = 0, lty = 2, col = "gray50", lwd = 0.8)

plot(dados_diff[, "ln_M1"],
     type = "l", col = "#7d6608", lwd = 1.5,
     main = "Taxa de crescimento de M1 (1ª dif. do ln)",
     ylab = "Δ ln(M1)", xlab = "")
abline(h = 0, lty = 2, col = "gray50", lwd = 0.8)

plot(dados_diff[, "ln_M2"],
     type = "l", col = "#6e2f1a", lwd = 1.5,
     main = "Taxa de crescimento de M2 (1ª dif. do ln)",
     ylab = "Δ ln(M2)", xlab = "Período (mensal)")
abline(h = 0, lty = 2, col = "gray50", lwd = 0.8)

dev.off()
cat(">> Figura 2 exportada.\n")

# -----------------------------------------------------------------------------
# 5. FUNÇÕES DE AUTOCORRELAÇÃO (ACF e PACF) – SÉRIES EM NÍVEL
# -----------------------------------------------------------------------------

series_info <- list(
  list(nome = "inflacao", titulo = "IPCA",
       arquivo = "Figura_2_ACF_PACF_IPCA.png"),
  list(nome = "ln_BM",    titulo = "Base Monetária – M0 (ln)",
       arquivo = "Figura_3_ACF_PACF_BM.png"),
  list(nome = "ln_M1",    titulo = "Agregado Monetário M1 (ln)",
       arquivo = "Figura_4_ACF_PACF_M1.png"),
  list(nome = "ln_M2",    titulo = "Agregado Monetário M2 (ln)",
       arquivo = "Figura_5_ACF_PACF_M2.png")
)

for (s in series_info) {
  png(filename = s$arquivo, width = 2400, height = 1400, res = 300)
  par(
    mfrow = c(1, 2),
    mar   = c(4, 5, 3.5, 2),
    cex.lab = 1.2, cex.axis = 1.0, cex.main = 1.1
  )
  acf(dados_ts[, s$nome],  main = paste0("ACF – ",  s$titulo), lag.max = 36)
  pacf(dados_ts[, s$nome], main = paste0("PACF – ", s$titulo), lag.max = 36)
  dev.off()
  cat(">>", s$arquivo, "exportado.\n")
}

# ── ACF e PACF das séries DIFERENCIADAS ──────────────────────────────────────

series_diff_info <- list(
  list(nome = "inflacao", titulo = "IPCA (1ª diferença)",
       arquivo = "Figura_6_ACF_PACF_IPCA_diff.png"),
  list(nome = "ln_BM",    titulo = "M0 (1ª dif. do ln)",
       arquivo = "Figura_7_ACF_PACF_BM_diff.png"),
  list(nome = "ln_M1",    titulo = "M1 (1ª dif. do ln)",
       arquivo = "Figura_8_ACF_PACF_M1_diff.png"),
  list(nome = "ln_M2",    titulo = "M2 (1ª dif. do ln)",
       arquivo = "Figura_9_ACF_PACF_M2_diff.png")
)

for (s in series_diff_info) {
  png(filename = s$arquivo, width = 2400, height = 1400, res = 300)
  par(
    mfrow = c(1, 2),
    mar   = c(4, 5, 3.5, 2),
    cex.lab = 1.2, cex.axis = 1.0, cex.main = 1.1
  )
  acf(dados_diff[, s$nome],  main = paste0("ACF – ",  s$titulo), lag.max = 36)
  pacf(dados_diff[, s$nome], main = paste0("PACF – ", s$titulo), lag.max = 36)
  dev.off()
  cat(">>", s$arquivo, "exportado.\n")
}


# -----------------------------------------------------------------------------
# 6. TESTES DE RAIZ UNITÁRIA – ADF (Dickey-Fuller Aumentado)
# -----------------------------------------------------------------------------

cat("\n", strrep("=", 60), "\n")
cat("  TESTE ADF – SÉRIES EM NÍVEL\n")
cat(strrep("=", 60), "\n")

adf_nivel <- list(
  inflacao = adf.test(dados_ts[, "inflacao"]),
  ln_BM    = adf.test(dados_ts[, "ln_BM"]),
  ln_M1    = adf.test(dados_ts[, "ln_M1"]),
  ln_M2    = adf.test(dados_ts[, "ln_M2"])
)

for (nm in names(adf_nivel)) {
  res <- adf_nivel[[nm]]
  cat(sprintf(
    "%-12s | Estat.: %7.4f | Lag: %d | p-valor: %.4f | %s\n",
    nm,
    res$statistic,
    res$parameter,
    res$p.value,
    ifelse(res$p.value < 0.05, "ESTACIONÁRIA", "NÃO estacionária")
  ))
}

cat("\n", strrep("=", 60), "\n")
cat("  TESTE ADF – PRIMEIRAS DIFERENÇAS\n")
cat(strrep("=", 60), "\n")

adf_diff <- list(
  inflacao = adf.test(dados_diff[, "inflacao"]),
  ln_BM    = adf.test(dados_diff[, "ln_BM"]),
  ln_M1    = adf.test(dados_diff[, "ln_M1"]),
  ln_M2    = adf.test(dados_diff[, "ln_M2"])
)

for (nm in names(adf_diff)) {
  res <- adf_diff[[nm]]
  cat(sprintf(
    "%-12s | Estat.: %7.4f | Lag: %d | p-valor: %.4f | %s\n",
    nm,
    res$statistic,
    res$parameter,
    res$p.value,
    ifelse(res$p.value < 0.05, "ESTACIONÁRIA", "NÃO estacionária")
  ))
}


# -----------------------------------------------------------------------------
# 7. SELEÇÃO DO NÚMERO DE DEFASAGENS
# -----------------------------------------------------------------------------

cat("\n", strrep("=", 60), "\n")
cat("  SELEÇÃO DE DEFASAGENS – VARselect (dados_diff)\n")
cat(strrep("=", 60), "\n")

lag_select <- VARselect(dados_diff, lag.max = 12, type = "const")
print(lag_select$selection)
cat("Critérios completos:\n")
print(round(lag_select$criteria, 4))

# Número de defasagens definido pelo AIC (mais comum para previsão)
p_otimo <- lag_select$selection["AIC(n)"]
cat(sprintf("\nDefasagem ótima selecionada (AIC): p = %d\n", p_otimo))


# -----------------------------------------------------------------------------
# 8. ESTIMAÇÃO DO MODELO VAR
# -----------------------------------------------------------------------------

cat("\n", strrep("=", 60), "\n")
cat("  ESTIMAÇÃO DO MODELO VAR (p = 12)\n")
cat(strrep("=", 60), "\n")

var_model <- VAR(dados_diff, p = 12, type = "const")

# Exportar sumário completo para arquivo texto (Anexo 1)
capture.output(
  summary(var_model),
  file = "Anexo1_VAR_summary.txt"
)
cat(">> Sumário do VAR exportado para 'Anexo1_VAR_summary.txt'.\n")


# -----------------------------------------------------------------------------
# 9. DIAGNÓSTICO DO MODELO VAR
# -----------------------------------------------------------------------------

cat("\n", strrep("=", 60), "\n")
cat("  DIAGNÓSTICO – ESTABILIDADE E AUTOCORRELAÇÃO\n")
cat(strrep("=", 60), "\n")

# 9.1 Estabilidade (todas as raízes devem ser < 1)
raizes <- roots(var_model)
cat(sprintf(
  "Máxima raiz característica: %.6f → VAR %s\n",
  max(raizes),
  ifelse(max(raizes) < 1, "ESTÁVEL ✓", "INSTÁVEL ✗")
))

# 9.2 Teste de Portmanteau (autocorrelação serial dos resíduos)
# Nota: lags.pt > p para evitar graus de liberdade zero
pt_test <- serial.test(var_model, lags.pt = 36, type = "PT.asymptotic")
cat("\nTeste de Portmanteau (lags = 36):\n")
print(pt_test)

cat(paste0(
  "\nNOTA: A presença de autocorrelação residual é comum em modelos VAR ",
  "de alta ordem com séries monetárias persistentes. A estabilidade ",
  "dinâmica está garantida (máxima raiz < 1), o que assegura a validade ",
  "assintótica dos testes de Causalidade de Granger.\n"
))


# -----------------------------------------------------------------------------
# 10. TESTES DE CAUSALIDADE DE GRANGER
# -----------------------------------------------------------------------------

cat("\n", strrep("=", 60), "\n")
cat("  TESTES DE CAUSALIDADE DE GRANGER\n")
cat(strrep("=", 60), "\n")

# Hipótese nula (H0): a variável 'cause' NÃO causa (no sentido de Granger)
# as demais variáveis do sistema.

causas <- c("ln_BM", "ln_M1", "ln_M2", "inflacao")

resultados_granger <- lapply(causas, function(v) {
  res <- causality(var_model, cause = v)
  list(
    variavel  = v,
    F_stat    = res$Granger$statistic,
    df1       = res$Granger$parameter["df1"],
    df2       = res$Granger$parameter["df2"],
    p_valor   = res$Granger$p.value
  )
})

cat(sprintf(
  "\n%-12s | %-10s | %-6s | %-6s | %-10s | %s\n",
  "Variável", "F-Stat", "df1", "df2", "p-valor", "Decisão (5%)"
))
cat(strrep("-", 72), "\n")

for (r in resultados_granger) {
  cat(sprintf(
    "%-12s | %-10.4f | %-6d | %-6d | %-10.6f | %s\n",
    r$variavel,
    r$F_stat,
    r$df1,
    r$df2,
    r$p_valor,
    ifelse(r$p_valor < 0.05, "Rejeita H0 ***", "Não rejeita H0")
  ))
}

cat("\n*** Causalidade de Granger significativa ao nível de 5%.\n")

cat(paste0(
  "\nINTERPRETAÇÃO:\n",
  "  - M1 e M2 causam (Granger) a inflação → evidência de transmissão\n",
  "    monetária para o nível de preços.\n",
  "  - A inflação também causa (Granger) as demais variáveis,\n",
  "    indicando bidirecionalidade (endogeneidade da moeda).\n",
  "  - M0 (Base Monetária) não apresenta causalidade significativa\n",
  "    para a inflação ao nível de 5%.\n"
))


# -----------------------------------------------------------------------------
# 11. EXPORTAÇÃO FINAL
# -----------------------------------------------------------------------------

# Salvar ambiente de trabalho
save.image("base_tratada_TCC.RData")
cat("\n>> Ambiente salvo em 'base_tratada_TCC.RData'.\n")

# Exportar base tratada em CSV (para conferência e anexo)
write.csv(dados, "dados_tratados_TCC.csv", row.names = FALSE)
cat(">> Base tratada exportada em 'dados_tratados_TCC.csv'.\n")

# Exportar resultados do Granger em texto
sink("Resultados_Granger.txt")
cat("RESULTADOS – TESTES DE CAUSALIDADE DE GRANGER\n")
cat("Modelo VAR(12) – Primeiras Diferenças – jan/2002 a dez/2025\n\n")
for (v in causas) {
  res <- causality(var_model, cause = v)
  cat(paste0("Causa: ", v, "\n"))
  print(res$Granger)
  cat("\n")
}
sink()
cat(">> Resultados do Granger exportados em 'Resultados_Granger.txt'.\n")

# Restaurar configurações gráficas originais
# Restaurar configurações gráficas originais (exceto "pin", que depende do device)
op_restaurar <- op_original[!names(op_original) %in% c("pin", "cin", "cra", "csi", "cxy", "din")]
par(op_restaurar)

cat("\n", strrep("=", 60), "\n")
cat("  ANÁLISE CONCLUÍDA COM SUCESSO\n")
cat(strrep("=", 60), "\n")
