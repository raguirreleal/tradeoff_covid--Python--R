# Tradeoff entre morte e desemprego na pandemia de COVID-19

Estudo da taxa de troca no tradeoff enfrentado, no início da pandemia de covid19, entre mortes e desemprego. 
As decisões governamentais pelo isolamento físico tinham que escolher entre ter:
- alto nível de isolamento com:
  - baixa quantidade de mortes; e
  - alta quantidade de desemprego.
- baixo nível de isolamento com:
  - alta quantidade de mortes; e
  - baixa quantidade de desemprego
A partir de um modelo epidemiológico Susceptible-Exposed-Infectious-Recovered (SEIR), em python, estimou-se a taxa de mortalidade para os diferentes níveis de isolamento, para todos os estados brasileiros, suas grandes regiões, e para o país.
Com um modelo de Vetores Auto Regressivos (VAR), em R, buscou-se estimar a relação entre a taxa de isolamento e a taxa de desemprego, para as mesmas localidades.
Em R foi feita a análise entre essas duas taxas.

***

Contém web scraping de dados de isolamento/mobilidade (repositório do Google), de vacinação, infecções, mortes etc. (repositório do DataSUS) e de desemprego (IBGE/PNAD).
Faz download e tratamento dos dados básicos e os salva para utilização no modelo SEIR.
Realiza simulação com o modelo SEIR, estimando os melhores parâmetros com otimização numérica (força bruta).
Faz análises dos dados com modelos ARIMA, dessazonalização com o X13-ARIMA/SEATS. Utiliza métodos não paramétricos de splines para desenhar curvas (funções) para as séries de tempo. Cria gráficos e salva tabelas para análise.

