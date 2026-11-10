using ITensors
using ITensorMPS

using Printf

using Random
Random.seed!()

# --- Reporte de Errores ---
struct QASMSyntaxError <: Exception
    mensaje::String
    linea::UInt64
end

Base.show(io::IO, e::QASMSyntaxError) = print(io, "Error en línea $(e.linea):\n\t$(e.mensaje)")

# --- Definiciones Básicas ---
const Qubit = UInt8

# Set de Instrucciones de la Arquitectura (ISA)
abstract type Instruccion end

# Instrucciones Simples
struct GateSimple <: Instruccion
    target::Qubit
    nombre::String
end

# Instrucciones de Control
struct GateControl <: Instruccion
    control::Qubit
    nombre::String
    target::Qubit
end

struct GateControl2 <: Instruccion
    control1::Qubit
    control2::Qubit
    target::Qubit
end

# Instrucciones de Rotación
struct GateRotation2 <: Instruccion
    target::Qubit
    nombre::String
    theta::Float64
end

struct GateRotation3 <: Instruccion
    target::Qubit
    theta::Float64
    phi::Float64
    lambda::Float64
end

# Instrucción de Medición
struct OpMedir <: Instruccion
    q_target::Qubit
    c_target::UInt8
end

# Instrucción de Impresión y DEBUG
struct OpPrint <: Instruccion
    c_target::UInt8
end

struct OpDebug <: Instruccion end

# Instrucción de Control de Flujo
struct OpJmp <: Instruccion
    line::UInt64
end

# Instrucciones de Control de Flujo Condicionado
 @enum Condicion IF_ZERO=0 IF_ONE=1 IF_EQ=2 IF_NEQ=3

struct OpJmpCond <: Instruccion
    control::UInt8
    line::UInt64
    condition::Condicion
end

struct OpJmpCond2 <: Instruccion
    control1::UInt8
    control2::UInt8
    line::UInt64
    condition::Condicion
end

# Instrucciones Soportadas
@enum TipoInstruccion X Y Z H S T RX RY RZ U CNOT CZ SWAP CCX MEDIR BARRERA GETQ GETC RESET PRINT DEBUG JMP JO JZ JE JD

# --- Definiciones para el Lexer ---
@enum TokenType Instruction QuBit ClassicBit Number Comma Newline Label EOF_TK

struct Token
    type::TokenType
    lexeme::String
    linea::UInt64
    value::Float64

end
# Constructor por defecto para value (0.0)
Token(type, lexeme, linea) = Token(type, lexeme, linea, 0.0)

# --- Estado Cuántico ---
mutable struct QuantumState
    indices::Vector{Index{Int}}
    psi::MPS

    QuantumState() = new([], MPS())
end

function init!(qs::QuantumState, total_qubits::Integer)
    # 1. Creamos los índices físicos (dimensión 2 para qubits)
    # "Site, Qubit" son etiquetas para poder buscarlos después
    qs.indices = siteinds("Qubit", total_qubits)

    # 2. Inicializamos el estado |000...0>
    # Esto crea un MPS de norma 1
    qs.psi = MPS(qs.indices, "0")
end

# --- Lexer (S -> T) ---
mutable struct Lexer
    source::String
    cursor::Int
    len::Int
    linea::UInt64

    Lexer(src::String) = new(src, 1, lastindex(src), 1)
end

function peek_char(l::Lexer) # Leer un caracter, si estamos al final del archivo, devolvemos el terminador nulo
    l.cursor > l.len ? '\0' : l.source[l.cursor]
end

function get_char!(l::Lexer) # Consumir un caracter y retornarlo
    c = peek_char(l)
    l.cursor = nextind(l.source, l.cursor)
    return c
end

function tokenize(l::Lexer) # Tokenizar
    tokens = Token[]
    while l.cursor <= l.len
        c = peek_char(l)

        if c in [' ', '\t', '\r'] # Ignoramos los espacios, tabulaciones y retornos de carro
            get_char!(l)

        elseif c == '@' # Etiquetas
            start = l.cursor
            get_char!(l) # Consumimos el "@"
            while isdigit(peek_char(l)) || isletter(peek_char(l)) get_char!(l) end
            push!(tokens, Token(Label, l.source[start:l.cursor - 1], l.linea))

        elseif c == '\n' # Saltos de línea
            push!(tokens, Token(Newline, "\n", l.linea))
            l.linea += 1
            get_char!(l)

        elseif c == ';' # Comentarios
            while peek_char(l) != '\n' && peek_char(l) != '\0' get_char!(l) end # Hasta no encontrar un salto de línea o "\0"

        elseif c == ',' # Comas
            push!(tokens, Token(Comma, ",", l.linea))
            get_char!(l)

        elseif isdigit(c) # Números
            push!(tokens, read_number!(l))

        elseif isletter(c) # Caracteres alfabéticos
            push!(tokens, read_variable!(l))

        else # Desconocido, se descarta
            get_char!(l)

        end
    end
    push!(tokens, Token(EOF_TK, "EOF", l.linea)) # Final del archivo, agregamos el token EOF (End Of File)
    return tokens
end

function read_number!(l::Lexer) # Función para leer números
    start = l.cursor
    while isdigit(peek_char(l)) || peek_char(l) == '.' get_char!(l) end
    lexeme = l.source[start:l.cursor - 1]
    val = parse(Float64, lexeme) # Se lee el número y almacena como un float

    if peek_char(l) == 'd' # El sufijo "d" al final de un número indica que se trata de Grados
        val = val * pi / 180 # Conversión a radianes
        get_char!(l) # Se consume la "d"
    end
    return Token(Number, lexeme, l.linea, val)
end

function read_variable!(l::Lexer) # Función para leer variables
    start = l.cursor
    while isdigit(peek_char(l)) || isletter(peek_char(l)) get_char!(l) end
    lexeme = l.source[start:l.cursor - 1] # Lee toda la cadena de caracteres

    # Caso 1 : Se trata de un Qubit (prefijo Q)
    if startswith(lexeme, "q") && length(lexeme) > 1 && all(isdigit, lexeme[2:end])
        return Token(QuBit, lexeme, l.linea, parse(Float64, lexeme[2:end]))
    # Caso 2 : Se trata de un Bit clásico (Prefijo c)
    elseif startswith(lexeme, "c") && length(lexeme) > 1 && all(isdigit, lexeme[2:end])
        return Token(ClassicBit, lexeme, l.linea, parse(Float64, lexeme[2:end]))
    end
    return Token(Instruction, lexeme, l.linea)
end

# --- Parser (T -> IR)---
mutable struct Program
    code::Vector{Instruccion}
    max_qubits_used::Int
    mapaEtiquetas::Dict{String, Int}

    Program() = new([], 0, Dict())
end

function validar_token!(tk::Token, tipo_esperado::TokenType, msg::String)
    if tk.type != tipo_esperado
        if tk.lexeme == "\n"
            throw(QASMSyntaxError("$msg (se encontró '\\n')", tk.linea))
        else
            throw(QASMSyntaxError("$msg (se encontró '$(tk.lexeme)')", tk.linea))
        end
    end
end

function parse_program(tokens::Vector{Token})
    prog = Program()
    pos = 1

    peek_tk() = pos <= length(tokens) ? tokens[pos] : Token(EOF_TK, "")
    get_tk!() = (t = peek_tk(); pos += 1; t)

    # Primera pasada: Guardamos las Etiquetas
    pos = 1
    dummy_counter = 1
    while pos <= length(tokens)
        t = tokens[pos]

        if t.type == Label
            # Registramos si la etiqueta no es un argumento
            prog.mapaEtiquetas[t.lexeme] = dummy_counter
            pos += 1
        elseif t.type == Instruction
            # Si es instrucción, incrementamos el contador de código real
            if !(t.lexeme in ["GETQ", "GETC"])
                dummy_counter += 1
            end
            pos += 1
            # Salteamos todos los tokens siguientes hasta
            # encontrar un Newline para NO procesar las etiquetas usadas como
            # arugmentos (Ej, JMP @etiqueta)
            while pos <= length(tokens) && tokens[pos].type != Newline && tokens[pos].type != EOF_TK
                pos += 1
            end
        else
            pos += 1
        end
    end

    # Segunda pasada: Generación de código
    pos = 1
    while peek_tk().type != EOF_TK # Hasta el final del archivo
        tk = get_tk!()
        if tk.type == Newline continue end
        if tk.type != Instruction continue end

        instr_str = tk.lexeme

        # Caso 1 : Puertas simples
        if instr_str in ["X", "Y", "Z", "H", "S", "T"]
            q = get_tk!()
            validar_token!(q, QuBit, "Se esperaba un identificador de Qubit")

            push!(prog.code, GateSimple(UInt8(q.value), instr_str))
            prog.max_qubits_used = max(prog.max_qubits_used, Int(q.value) + 1)

        # Caso 2 : Instrucción de Medición
        elseif instr_str == "MEDIR"
            q = get_tk!()
            validar_token!(q, QuBit, "Se esperaba un identificador de Qubit")
            comma_tk = get_tk!()
            validar_token!(comma_tk, Comma, "Se esperaba una coma")
            c = get_tk!()
            validar_token!(c, ClassicBit, "Se esperaba un identificador de Bit")

            push!(prog.code, OpMedir(UInt8(q.value), UInt8(c.value)))
            prog.max_qubits_used = max(prog.max_qubits_used, Int(q.value) + 1)

        # Caso 3 : Control Simple e Intercambios
        elseif instr_str in ["CNOT", "CZ", "SWAP"]
            c = get_tk!() # Qubit de control
            validar_token!(c, QuBit, "Se esperaba un identificador de Qubit")
            comma_tk = get_tk!()
            validar_token!(comma_tk, Comma, "Se esperaba una coma")
            t = get_tk!()
            validar_token!(t, QuBit, "Se esperaba un identificador de Qubit")

            push!(prog.code, GateControl(UInt8(c.value), instr_str, UInt8(t.value)))
            prog.max_qubits_used = max(prog.max_qubits_used, Int(c.value) + 1, Int(t.value) + 1)

        # Caso 4 : Control Doble (Toffoli)
        elseif instr_str == "CCX"
            q1 = get_tk!()
            validar_token!(q1, QuBit, "Se esperaba un identificador de Qubit")
            comma_tk = get_tk!()
            validar_token!(comma_tk, Comma, "Se esperaba una coma")
            q2 = get_tk!()
            validar_token!(q2, QuBit, "Se esperaba un identificador de Qubit")
            comma_tk = get_tk!()
            validar_token!(comma_tk, Comma, "Se esperaba una coma")
            t = get_tk!()
            validar_token!(t, QuBit, "Se esperaba un identificador de Qubit")

            push!(prog.code, GateControl2(UInt8(q1.value), UInt8(q2.value), UInt8(t.value)))
            prog.max_qubits_used = max(prog.max_qubits_used, Int(q1.value) + 1, Int(q2.value) + 1, Int(t.value) + 1)

        # -- Rotaciones  --
        # Caso 5 : Rotación Simple
        elseif instr_str in ["RX", "RY", "RZ"]
            q = get_tk!()
            validar_token!(q, QuBit, "Se esperaba un identificador de Qubit")
            validar_token!(get_tk!(), Comma, "Se esperaba una coma")
            theta = get_tk!()
            validar_token!(theta, Number, "Se esperaba un número")
            name = titlecase(instr_str)
            push!(prog.code, GateRotation2(UInt8(q.value), name, theta.value))

        # Caso 6 : Rotación Unitaria
        elseif instr_str == "U"
            q = get_tk!()
            validar_token!(q, QuBit, "Se esperaba un identificador de Qubit")
            args = Float64[]
            for i in 1:3
                validar_token!(get_tk!(), Comma, "Se esperaba una coma antes del parámetro $i")
                numero = get_tk!()
                validar_token!(numero, Number, "Se esperaba un valor numérico para el parámetro $i")
                push!(args, numero.value)
            end
            push!(prog.code, GateRotation3(UInt8(q.value), args[1], args[2], args[3]))
            prog.max_qubits_used = max(prog.max_qubits_used, Int(q.value) + 1)

        # -- Saltos --
        # Caso 7 : Salto Incondicional
        elseif instr_str == "JMP"
            target = get_tk!()
            if !haskey(prog.mapaEtiquetas, target.lexeme)
                throw(QASMSyntaxError("La etiqueta '$(target.lexeme)' no existe", target.linea))
            end
            push!(prog.code, OpJmp(prog.mapaEtiquetas[target.lexeme]))

        # Caso 8 : Salto Condicional Simple (Jump If Zero y Jump If One)
        elseif instr_str in ["JZ", "JO"]
            cond_bit = get_tk!()
            validar_token!(cond_bit, ClassicBit, "Se esperaba un identificador de Bit")
            validar_token!(get_tk!(), Comma, "Se esperaba una coma")
            label = get_tk!()

            if !haskey(prog.mapaEtiquetas, label.lexeme)
                throw(QASMSyntaxError("La etiqueta '$(label.lexeme)' no existe", label.linea))
            end
            cond_val = (instr_str == "JZ") ? IF_ZERO : IF_ONE
            push!(prog.code, OpJmpCond(UInt8(cond_bit.value), prog.mapaEtiquetas[label.lexeme], cond_val))

        # Caso 9 : Salto Condicional Compuesto (Jump If Equal y Jump If Not Equal)
        elseif instr_str in ["JE", "JD"]
            c1 = get_tk!()
            validar_token!(c1, ClassicBit, "Se esperaba un identificador de Bit")
            validar_token!(get_tk!(), Comma, "Se esperaba una coma")
            c2 = get_tk!()
            validar_token!(c2, ClassicBit, "Se esperaba un identificador de Bit")
            validar_token!(get_tk!(), Comma, "Se esperaba una coma")
            label = get_tk!()

             if !haskey(prog.mapaEtiquetas, label.lexeme)
                throw(QASMSyntaxError("La etiqueta '$(label.lexeme)' no existe", label.linea))
            end
            cond_val = (instr_str == "JE") ? IF_EQ : IF_NEQ
            push!(prog.code, OpJmpCond2(UInt8(c1.value), UInt8(c2.value), prog.mapaEtiquetas[label.lexeme], cond_val))

        # Caso 10 : Mostrar un Bit clásico
        elseif instr_str == "PRINT"
            c = get_tk!()
            push!(prog.code, OpPrint(UInt8(c.value)))

        # Caso 11 : Información del Sistema
        elseif instr_str == "DEBUG"
            push!(prog.code, OpDebug())

        end
    end
    return prog
end

# --- Máquina Virtual (IR -> Ejecución) ---
function ejecutar_vm(prog::Program, cutoff_val::Float64, maxdim_val::Int, topology_limit::Int)
    println("Inicializando sistema de $(prog.max_qubits_used) qubit(s)...")

    # Inicialización ITensors
    s = siteinds("Qubit", prog.max_qubits_used)
    psi = MPS(s, "0")

    # Memoria clásica
    classical_memory = Dict{UInt8, Bool}()

    println("Ejecutando...")

    ip = 1
    while ip <= length(prog.code)
        inst = prog.code[ip]
        target_idx = 0

        if inst isa GateSimple
            g = op(inst.nombre, s[inst.target + 1])
            psi = apply(g, psi; cutoff=cutoff_val, maxdim=maxdim_val)

        elseif inst isa OpMedir
            # 1. Preparación
            orthogonalize!(psi, inst.q_target + 1)

            # 2. Calcular la probabilidad de obtener ∣1⟩
            p1_op = op("Proj1", s[inst.q_target + 1])
            prob1 = expect(psi, "Proj1"; sites=inst.q_target + 1)

            # 3. Obtenemos un resultado
            # Número pseudo aleatorio entre 0 y 1. Si es menor que la prob de ∣1⟩, el resultado es 1.
            res = (rand() < prob1) ? 1 : 0

            # 4. Guardar en la memoria clásica (Boolean)
            classical_memory[inst.c_target] = (res == 1)

            # 5. Aplicar el colapso (Proyectar el estado)
            # Si salió 0, aplicamos Proj0; si salió 1, Proj1.
            nombre_proyector = (res == 0) ? "Proj0" : "Proj1"
            P = op(nombre_proyector, s[inst.q_target + 1])

            # Aplicamos el proyector y normalizamos
            psi = apply(P, psi; cutoff=cutoff_val, maxdim=maxdim_val)
            normalize!(psi)

            orthogonalize!(psi, inst.q_target + 1)

        # --- Control ---
        elseif inst isa GateControl
            g = op(inst.nombre, s[inst.control + 1], s[inst.target + 1])
            psi = apply(g, psi; cutoff=cutoff_val, maxdim=maxdim_val)

        # --- Puerta de Control de 3 Qubits (CCX / Toffoli) ---
        elseif inst isa GateControl2
            g = op("CCNOT", s[inst.control1 + 1], s[inst.control2 + 1], s[inst.target + 1])
            psi = apply(g, psi; cutoff=cutoff_val, maxdim=maxdim_val)

        # --- Rotación ---
        elseif inst isa GateRotation2
            g = op(inst.nombre, s[inst.target + 1]; θ=inst.theta)
            psi = apply(g, psi; cutoff=cutoff_val, maxdim=maxdim_val)

        elseif inst isa GateRotation3 # U / Unitary
            # Descomposición de Euler
            # Aplicamos la secuencia: Rz(λ) -> Ry(θ) -> Rz(ϕ)

            # 1. Rotación en Z con ángulo lambda
            psi = apply(op("Rz", s[inst.target + 1]; θ=inst.lambda), psi; cutoff=cutoff_val, maxdim=maxdim_val)

            # 2. Rotación en Y con ángulo theta
            psi = apply(op("Ry", s[inst.target + 1]; θ=inst.theta), psi; cutoff=cutoff_val, maxdim=maxdim_val)

            # 3. Rotación en Z con ángulo phi
            psi = apply(op("Rz", s[inst.target + 1]; θ=inst.phi), psi; cutoff=cutoff_val, maxdim=maxdim_val)

        # --- Saltos ---
        # Incondicionales
        elseif inst isa OpJmp
            ip = inst.line
            continue

        # Condicionales Simples
        elseif inst isa OpJmpCond
            val = get(classical_memory, inst.control, false)
            met_condition = false
            if inst.condition == IF_ZERO && val == false met_condition = true end
            if inst.condition == IF_ONE  && val == true  met_condition = true end

            if met_condition
                ip = inst.line
                continue
            end

        # Condicionales Compuestos
        elseif inst isa OpJmpCond2
            val1 = get(classical_memory, inst.control1, false)
            val2 = get(classical_memory, inst.control2, false)
            met_condition = false

            if inst.condition == IF_EQ  && val1 == val2 met_condition = true end
            if inst.condition == IF_NEQ && val1 != val2 met_condition = true end

            if met_condition
                ip = inst.line
                continue
            end

        # Información de Debug
        elseif inst isa OpDebug
                println("\n" * "="^30)
                println()
                println("DEBUG: Estado del Sistema")
                println("Qubits activos: $(prog.max_qubits_used)")

                println("\nProyecciones de Bloch por Qubit:")
                for i in 1:prog.max_qubits_used
                    orthogonalize!(psi, i)

                    # 1. Extraemos las probabilidades
                    p1 = expect(psi, "Proj1"; sites=i)

                    # 2. Extraemos los valores esperados de X e Y para calcular la fase
                    ex = expect(psi, "X"; sites=i)
                    ey = expect(psi, "Y"; sites=i)

                    # 3. Calculamos la fase relativa en grados
                    # atan2(y, x) nos da el ángulo en el plano complejo
                    fase_rad = atan(ey, ex)
                    fase_deg = rad2deg(fase_rad)

                    # 4. Barra visual
                    bar_len = 15
                    filled = round(Int, p1 * bar_len)
                    bar = "█"^filled * "░"^(bar_len - filled)

                    @printf("  q%d: [%s] P( ∣1⟩ ):%5.1f%% | Fase: %6.1f°\n",
                            i-1, bar, p1*100, fase_deg)

                end

                # ITensors nos permite imprimir el MPS directamente para ver sus dimensiones
                println("\nEstructura del MPS:")
                show(stdout, "text/plain", psi)

                # Imprimir la memoria clásica si no está vacía
                if !isempty(classical_memory)
                    println("\nMemoria Clásica:")
                    for k in sort(collect(keys(classical_memory)))
                        println("  c$k : $(classical_memory[k] ? 1 : 0)")
                    end
                end
                println()
                println("="^30 * "\n")

        # Impresión
        elseif inst isa OpPrint
            val = get(classical_memory, inst.c_target, false)
            println("c$(inst.c_target) = $(val ? 1 : 0)")
        end # elseif

        ip += 1
    end # while
end # function

# --- Función Principal ---
function main()

    if length(ARGS) < 1
        println("Uso: julia qasm.jl <archivo.fq> [opciones]")
        println("\nOpciones:")
        println("  -cutoff=N    Exponente del cutoff (ej: 12 para 1e-12)")
        println("  -maxdim=N    Dimensión máxima de enlace (bond dimension)")
        println("  -topology=N  Distancia máxima permitida para puertas de 2 qubits")
        return
    end

    # Valores por defecto
    archivo = ""
    cutoff = 0.0
    maxdim = 1000000
    topology = 0

    # Procesar argumentos
    for arg in ARGS
        if startswith(arg, "-cutoff=")
            val = parse(Float64, split(arg, "=")[2])
            cutoff = 10.0^(-val)
        elseif startswith(arg, "-maxdim=")
            maxdim = parse(Int, split(arg, "=")[2])
        elseif startswith(arg, "-topology=")
            topology = parse(Int, split(arg, "=")[2])
        elseif !startswith(arg, "-")
            archivo = arg
        end
    end

    if isempty(archivo)
        println("Error: No se especificó un archivo .fq")
        return
    end

    try
        contenido = read(ARGS[1], String)
        lexer = Lexer(contenido)
        tokens = tokenize(lexer)
        prog = parse_program(tokens)
        ejecutar_vm(prog, cutoff, maxdim, topology)
    catch e
        if e isa QASMSyntaxError # Si es un error de nuestro t ipo
            println(e)
        elseif e isa SystemError
            println("Error: No se pudo abrir el archivo '$archivo'")
        else # Si es otro error, que se lance normalmente
            rethrow(e)
        end
    end
end

main()
