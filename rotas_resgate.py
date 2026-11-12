import heapq
import math
import random
import time

TAMANHO_MAPA = 20

PONTOS_RESGATE = {
    "Base - Hospital Central": (10, 10),
    "Ponto 1 - Abrigo Norte": (3, 17),
    "Ponto 2 - Escola Municipal": (15, 18),
    "Ponto 3 - Ginásio Poliesportivo": (18, 11),
    "Ponto 4 - Posto de Saúde Sul": (16, 3),
    "Ponto 5 - Comunidade Ribeirinha": (9, 2),
    "Ponto 6 - Vila Operária": (2, 4),
    "Ponto 7 - Praça da Matriz": (3, 10),
    "Ponto 8 - Hospital de Campanha": (12, 16),
}

NOMES_PONTOS = list(PONTOS_RESGATE.keys())
INDICES_ENTREGAS = list(range(1, len(PONTOS_RESGATE)))

MAPA = [[0 for _ in range(TAMANHO_MAPA)] for _ in range(TAMANHO_MAPA)]

for x in range(11, 16):
    MAPA[14][x] = 1
for y in range(11, 15):
    MAPA[y][14] = 1

for x in range(4, 9):
    MAPA[6][x] = 1
for y in range(4, 8):
    MAPA[y][6] = 1

for y in range(8, 13):
    MAPA[y][7] = 1


def heuristica_octile(a, b):
    dx = abs(a[0] - b[0])
    dy = abs(a[1] - b[1])
    return (dx + dy) + (math.sqrt(2) - 2) * min(dx, dy)


def busca_a_estrela_caminho(inicio, destino, mapa):
    movimentos = [
        (0, 1, 1.0), (0, -1, 1.0), (1, 0, 1.0), (-1, 0, 1.0),
        (1, 1, math.sqrt(2)), (1, -1, math.sqrt(2)),
        (-1, 1, math.sqrt(2)), (-1, -1, math.sqrt(2)),
    ]

    fronteira = []
    heapq.heappush(fronteira, (0.0, inicio))

    veio_de = {inicio: None}
    custo_g = {inicio: 0.0}
    nos_expandidos = 0

    while fronteira:
        _, atual = heapq.heappop(fronteira)
        nos_expandidos += 1

        if atual == destino:
            caminho = []
            curr = destino
            while curr is not None:
                caminho.append(curr)
                curr = veio_de[curr]
            caminho.reverse()
            return custo_g[destino], caminho, nos_expandidos

        for dx, dy, custo_passo in movimentos:
            nx, ny = atual[0] + dx, atual[1] + dy
            if 0 <= nx < TAMANHO_MAPA and 0 <= ny < TAMANHO_MAPA:
                if mapa[ny][nx] == 1:
                    continue

                novo_custo = custo_g[atual] + custo_passo
                vizinho = (nx, ny)

                if vizinho not in custo_g or novo_custo < custo_g[vizinho]:
                    custo_g[vizinho] = novo_custo
                    prioridade_f = novo_custo + heuristica_octile(vizinho, destino)
                    heapq.heappush(fronteira, (prioridade_f, vizinho))
                    veio_de[vizinho] = atual

    return float("inf"), [], nos_expandidos


print("Calculando matriz de rotas seguras com A*...")
NUM_PONTOS = len(NOMES_PONTOS)
MATRIZ_DISTANCIAS_A_ESTRELA = [[0.0] * NUM_PONTOS for _ in range(NUM_PONTOS)]
TOTAL_NOS_A_ESTRELA_TRECHOS = 0

for i in range(NUM_PONTOS):
    for j in range(i + 1, NUM_PONTOS):
        p1 = PONTOS_RESGATE[NOMES_PONTOS[i]]
        p2 = PONTOS_RESGATE[NOMES_PONTOS[j]]
        custo, _, expandidos = busca_a_estrela_caminho(p1, p2, MAPA)
        MATRIZ_DISTANCIAS_A_ESTRELA[i][j] = custo
        MATRIZ_DISTANCIAS_A_ESTRELA[j][i] = custo
        TOTAL_NOS_A_ESTRELA_TRECHOS += expandidos


def calcular_custo_rota(rota_entregas):
    distancia = MATRIZ_DISTANCIAS_A_ESTRELA[0][rota_entregas[0]]
    for i in range(len(rota_entregas) - 1):
        distancia += MATRIZ_DISTANCIAS_A_ESTRELA[rota_entregas[i]][rota_entregas[i + 1]]
    distancia += MATRIZ_DISTANCIAS_A_ESTRELA[rota_entregas[-1]][0]
    return distancia


def resolver_sequencia_com_a_estrela():
    todos_pontos = set(INDICES_ENTREGAS)
    menor_aresta = min(
        MATRIZ_DISTANCIAS_A_ESTRELA[i][j]
        for i in range(NUM_PONTOS)
        for j in range(NUM_PONTOS)
        if i != j
    )

    fronteira = []
    h_inicial = len(todos_pontos) * menor_aresta
    heapq.heappush(fronteira, (h_inicial, 0.0, 0, frozenset(), []))

    melhor_g = {}
    estados_expandidos = 0

    while fronteira:
        f, g, atual, visitados, caminho = heapq.heappop(fronteira)
        estados_expandidos += 1

        if len(visitados) == len(todos_pontos):
            custo_final = g + MATRIZ_DISTANCIAS_A_ESTRELA[atual][0]
            return custo_final, caminho, estados_expandidos

        estado_chave = (atual, visitados)
        if estado_chave in melhor_g and melhor_g[estado_chave] <= g:
            continue
        melhor_g[estado_chave] = g

        nao_visitados = todos_pontos - visitados
        for prox in nao_visitados:
            novo_g = g + MATRIZ_DISTANCIAS_A_ESTRELA[atual][prox]
            novos_visitados = visitados | {prox}
            faltam = len(todos_pontos) - len(novos_visitados) + 1
            novo_h = faltam * menor_aresta
            heapq.heappush(
                fronteira,
                (novo_g + novo_h, novo_g, prox, novos_visitados, caminho + [prox]),
            )

    return float("inf"), [], estados_expandidos


TAMANHO_POPULACAO = 50
TAXA_CRUZAMENTO = 0.85
TAXA_MUTACAO = 0.15
NUM_GERACOES = 80


def selecao_torneio(populacao, k=3):
    participantes = random.sample(populacao, k)
    participantes.sort(key=calcular_custo_rota)
    return participantes[0]


def cruzamento_ordenado(pai1, pai2):
    if random.random() > TAXA_CRUZAMENTO:
        return pai1.copy()

    tamanho = len(pai1)
    ponto1, ponto2 = sorted(random.sample(range(tamanho), 2))

    filho = [None] * tamanho
    filho[ponto1 : ponto2 + 1] = pai1[ponto1 : ponto2 + 1]

    posicao_atual = (ponto2 + 1) % tamanho
    for gene in pai2[ponto2 + 1 :] + pai2[: ponto2 + 1]:
        if gene not in filho:
            filho[posicao_atual] = gene
            posicao_atual = (posicao_atual + 1) % tamanho

    return filho


def mutacao_troca(rota):
    rota_mutada = rota.copy()
    if random.random() < TAXA_MUTACAO:
        i, j = random.sample(range(len(rota_mutada)), 2)
        rota_mutada[i], rota_mutada[j] = rota_mutada[j], rota_mutada[i]
    return rota_mutada


def resolver_sequencia_com_ag(semente=42):
    random.seed(semente)
    populacao = [
        random.sample(INDICES_ENTREGAS, len(INDICES_ENTREGAS))
        for _ in range(TAMANHO_POPULACAO)
    ]

    melhor_rota_global = None
    menor_distancia_global = float("inf")
    historico_geracoes = []
    avaliacoes = 0

    for geracao in range(1, NUM_GERACOES + 1):
        populacao.sort(key=calcular_custo_rota)
        avaliacoes += len(populacao)
        melhor_da_geracao = populacao[0]
        distancia_atual = calcular_custo_rota(melhor_da_geracao)

        if distancia_atual < menor_distancia_global:
            menor_distancia_global = distancia_atual
            melhor_rota_global = melhor_da_geracao.copy()

        if geracao % 10 == 0 or geracao == 1:
            historico_geracoes.append((geracao, menor_distancia_global, melhor_rota_global.copy()))

        nova_populacao = [populacao[0].copy(), populacao[1].copy()]

        while len(nova_populacao) < TAMANHO_POPULACAO:
            pai1 = selecao_torneio(populacao)
            pai2 = selecao_torneio(populacao)
            filho = cruzamento_ordenado(pai1, pai2)
            filho = mutacao_troca(filho)
            nova_populacao.append(filho)

        populacao = nova_populacao

    return menor_distancia_global, melhor_rota_global, avaliacoes, historico_geracoes


def exibir_mapa_ascii():
    print("\\n=== MAPA DA ÁREA DE DESASTRE (20x20) ===")
    print("Legenda: [B] Base | [1-8] Pontos de Resgate | [##] Zona Interditada | [ .] Livre\\n")
    coords_invertidas = {v: str(k) if k > 0 else "B" for k, (_, v) in enumerate(PONTOS_RESGATE.items())}

    for y in range(TAMANHO_MAPA - 1, -1, -1):
        linha = f"{y:2d} | "
        for x in range(TAMANHO_MAPA):
            if (x, y) in coords_invertidas:
                linha += f" {coords_invertidas[(x, y)]} "
            elif MAPA[y][x] == 1:
                linha += "## "
            else:
                linha += " . "
        print(linha)
    print("   +" + "---" * TAMANHO_MAPA)
    print("     " + "".join(f"{x:2d} " for x in range(TAMANHO_MAPA)))


if __name__ == "__main__":
    exibir_mapa_ascii()

    t0 = time.perf_counter()
    custo_astar, rota_astar, estados_astar = resolver_sequencia_com_a_estrela()
    tempo_astar = (time.perf_counter() - t0) * 1000

    t0 = time.perf_counter()
    custo_ag, rota_ag, avaliacoes_ag, historico_ag = resolver_sequencia_com_ag()
    tempo_ag = (time.perf_counter() - t0) * 1000

    print("\\n=== EVOLUÇÃO DO ALGORITMO GENÉTICO ===")
    print(f"{'Geração':<10} | {'Melhor Distância':<18} | {'Cromossomo (Ordem)'}")
    print("-" * 65)
    for ger, dist, crom in historico_ag:
        print(f"{ger:<10} | {dist:<18.2f} | {crom}")

    print("\\n=== QUADRO COMPARATIVO: A* vs ALGORITMO GENÉTICO ===")
    print(f"{'Métrica':<28} | {'A* (Busca Exata)':<20} | {'Algoritmo Genético':<20}")
    print("-" * 74)
    print(f"{'Distância Total da Rota':<28} | {custo_astar:<20.2f} | {custo_ag:<20.2f}")
    print(f"{'Estados / Avaliações':<28} | {estados_astar:<20} | {avaliacoes_ag:<20}")
    print(f"{'Tempo de Execução (ms)':<28} | {tempo_astar:<20.2f} | {tempo_ag:<20.2f}")
    print(f"{'Ordem dos Índices':<28} | {str(rota_astar):<20} | {str(rota_ag):<20}")

    print("\\n=== ITINERÁRIO DETALHADO (MELHOR ROTA) ===")
    roteiro = ["Base - Hospital Central"] + [NOMES_PONTOS[i] for i in rota_astar] + ["Base - Hospital Central"]
    for idx, parada in enumerate(roteiro):
        print(f" {idx}. {parada}")
