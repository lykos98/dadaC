	.file	"kdtree.c"
# GNU C17 (Ubuntu 13.2.0-4ubuntu3) version 13.2.0 (x86_64-linux-gnu)
#	compiled by GNU C version 13.2.0, GMP version 6.3.0, MPFR version 4.2.1, MPC version 1.3.1, isl version isl-0.26-GMP

# GGC heuristics: --param ggc-min-expand=100 --param ggc-min-heapsize=131072
# options passed: -mavx2 -mtune=generic -march=x86-64 -ggdb -O3 -fasynchronous-unwind-tables -fstack-protector-strong -fstack-clash-protection -fcf-protection
	.text
.Ltext0:
	.file 0 "/home/francesco/Desktop/dssc/robavaria/cc" "src/kdtree.c"
	.p2align 4
	.globl	swap
	.type	swap, @function
swap:
.LVL0:
.LFB53:
	.file 1 "src/kdtree.c"
	.loc 1 11 22 view -0
	.cfi_startproc
	.loc 1 11 22 is_stmt 0 view .LVU1
	endbr64	
	.loc 1 12 5 is_stmt 1 view .LVU2
	.loc 1 13 5 view .LVU3
.LVL1:
.LBB60:
.LBI60:
	.file 2 "/usr/include/x86_64-linux-gnu/bits/string_fortified.h"
	.loc 2 26 1 view .LVU4
.LBB61:
	.loc 2 29 3 view .LVU5
# /usr/include/x86_64-linux-gnu/bits/string_fortified.h:29:   return __builtin___memcpy_chk (__dest, __src, __len,
	.loc 2 29 10 is_stmt 0 discriminator 1 view .LVU6
	movq	(%rdi), %rax	# MEM <unsigned long> [(char * {ref-all})a_2(D)], _3
.LVL2:
	.loc 2 29 10 discriminator 1 view .LVU7
.LBE61:
.LBE60:
	.loc 1 14 5 is_stmt 1 view .LVU8
.LBB62:
.LBI62:
	.loc 2 26 1 view .LVU9
.LBB63:
	.loc 2 29 3 view .LVU10
# /usr/include/x86_64-linux-gnu/bits/string_fortified.h:29:   return __builtin___memcpy_chk (__dest, __src, __len,
	.loc 2 29 10 is_stmt 0 discriminator 1 view .LVU11
	movq	(%rsi), %rdx	# MEM <unsigned long> [(char * {ref-all})b_4(D)], _11
	movq	%rdx, (%rdi)	# _11, MEM <unsigned long> [(char * {ref-all})a_2(D)]
.LVL3:
	.loc 2 29 10 discriminator 1 view .LVU12
.LBE63:
.LBE62:
	.loc 1 15 5 is_stmt 1 view .LVU13
.LBB64:
.LBI64:
	.loc 2 26 1 view .LVU14
.LBB65:
	.loc 2 29 3 view .LVU15
# /usr/include/x86_64-linux-gnu/bits/string_fortified.h:29:   return __builtin___memcpy_chk (__dest, __src, __len,
	.loc 2 29 10 is_stmt 0 discriminator 1 view .LVU16
	movq	%rax, (%rsi)	# _3, MEM <unsigned long> [(char * {ref-all})b_4(D)]
.LVL4:
	.loc 2 29 10 discriminator 1 view .LVU17
.LBE65:
.LBE64:
	.loc 1 16 5 is_stmt 1 view .LVU18
# src/kdtree.c:17: }
	.loc 1 17 1 is_stmt 0 view .LVU19
	ret	
	.cfi_endproc
.LFE53:
	.size	swap, .-swap
	.p2align 4
	.globl	euclidean_distance
	.type	euclidean_distance, @function
euclidean_distance:
.LVL5:
.LFB54:
	.loc 1 19 62 is_stmt 1 view -0
	.cfi_startproc
	.loc 1 19 62 is_stmt 0 view .LVU21
	endbr64	
	.loc 1 20 5 is_stmt 1 view .LVU22
.LVL6:
	.loc 1 21 5 view .LVU23
.LBB66:
	.loc 1 21 9 view .LVU24
	.loc 1 21 30 discriminator 1 view .LVU25
	movl	data_dims(%rip), %ecx	# data_dims, data_dims.0_19
.LBE66:
# src/kdtree.c:19: FLOAT_TYPE euclidean_distance(FLOAT_TYPE* p1, FLOAT_TYPE* p2){
	.loc 1 19 62 is_stmt 0 view .LVU26
	movq	%rdi, %rdx	# tmp150, p1
.LBB70:
# src/kdtree.c:21:     for(unsigned int i = 0; i<data_dims; ++i){
	.loc 1 21 30 discriminator 1 view .LVU27
	testl	%ecx, %ecx	# data_dims.0_19
	je	.L11	#,
	leal	-1(%rcx), %eax	#, tmp128
	cmpl	$2, %eax	#, tmp128
	jbe	.L12	#,
	movl	%ecx, %edi	# data_dims.0_19, bnd.29
.LVL7:
	.loc 1 21 30 discriminator 1 view .LVU28
	xorl	%eax, %eax	# ivtmp.60
.LBE70:
# src/kdtree.c:20:     FLOAT_TYPE d = 0;
	.loc 1 20 16 view .LVU29
	vxorpd	%xmm0, %xmm0, %xmm0	# <retval>
	shrl	$2, %edi	#,
	salq	$5, %rdi	#, _87
.LVL8:
	.p2align 4,,10
	.p2align 3
.L6:
.LBB71:
.LBB67:
	.loc 1 22 9 is_stmt 1 view .LVU30
# src/kdtree.c:22:         FLOAT_TYPE dd = (p1[i] - p2[i]);
	.loc 1 22 20 is_stmt 0 view .LVU31
	vmovupd	(%rdx,%rax), %ymm4	# MEM <vector(4) double> [(double *)p1_13(D) + ivtmp.60_75 * 1], tmp154
	vsubpd	(%rsi,%rax), %ymm4, %ymm1	# MEM <vector(4) double> [(double *)p2_14(D) + ivtmp.60_75 * 1], tmp154, vect_dd_15.38
.LVL9:
	.loc 1 23 9 is_stmt 1 view .LVU32
	addq	$32, %rax	#, ivtmp.60
# src/kdtree.c:23:         d += dd*dd;
	.loc 1 23 16 is_stmt 0 view .LVU33
	vmulpd	%ymm1, %ymm1, %ymm1	# vect_dd_15.38, vect_dd_15.38, vect__8.39
	vaddsd	%xmm1, %xmm0, %xmm0	# stmp_d_16.40, <retval>, stmp_d_16.40
.LVL10:
	.loc 1 23 16 view .LVU34
	vunpckhpd	%xmm1, %xmm1, %xmm2	# tmp133, stmp_d_16.40
	vextractf128	$0x1, %ymm1, %xmm1	# vect__8.39, tmp135
	vaddsd	%xmm2, %xmm0, %xmm0	# stmp_d_16.40, stmp_d_16.40, stmp_d_16.40
# src/kdtree.c:23:         d += dd*dd;
	.loc 1 23 11 view .LVU35
	vaddsd	%xmm1, %xmm0, %xmm0	# stmp_d_16.40, stmp_d_16.40, stmp_d_16.40
	vunpckhpd	%xmm1, %xmm1, %xmm1	# tmp135, stmp_d_16.40
	vaddsd	%xmm1, %xmm0, %xmm0	# stmp_d_16.40, stmp_d_16.40, <retval>
.LVL11:
	.loc 1 23 11 view .LVU36
.LBE67:
	.loc 1 21 42 is_stmt 1 discriminator 3 view .LVU37
	.loc 1 21 30 discriminator 1 view .LVU38
	cmpq	%rdi, %rax	# _87, ivtmp.60
	jne	.L6	#,
	testb	$3, %cl	#, data_dims.0_19
	je	.L23	#,
	movl	%ecx, %eax	# data_dims.0_19, tmp.44
	andl	$-4, %eax	#,
	vzeroupper
.LVL12:
.L5:
	.loc 1 21 30 is_stmt 0 discriminator 1 view .LVU39
	subl	%eax, %ecx	# tmp.44, niters.41
	cmpl	$1, %ecx	#, niters.41
	je	.L9	#,
	movl	%eax, %edi	# tmp.44, tmp.44
.LVL13:
.LBB68:
	.loc 1 22 9 is_stmt 1 view .LVU40
# src/kdtree.c:22:         FLOAT_TYPE dd = (p1[i] - p2[i]);
	.loc 1 22 20 is_stmt 0 view .LVU41
	vmovupd	(%rdx,%rdi,8), %xmm5	# MEM <vector(2) double> [(double *)vectp_p1.46_83], tmp156
	vsubpd	(%rsi,%rdi,8), %xmm5, %xmm1	# MEM <vector(2) double> [(double *)vectp_p2.49_89], tmp156, vect_dd_27.51
.LVL14:
	.loc 1 23 9 is_stmt 1 view .LVU42
# src/kdtree.c:23:         d += dd*dd;
	.loc 1 23 16 is_stmt 0 view .LVU43
	vmulpd	%xmm1, %xmm1, %xmm1	# vect_dd_27.51, vect_dd_27.51, vect__26.52
	vaddsd	%xmm1, %xmm0, %xmm0	# stmp_d_25.53, <retval>, stmp_d_25.53
.LVL15:
# src/kdtree.c:23:         d += dd*dd;
	.loc 1 23 11 view .LVU44
	vunpckhpd	%xmm1, %xmm1, %xmm1	# vect__26.52, stmp_d_25.53
	vaddsd	%xmm0, %xmm1, %xmm0	# stmp_d_25.53, stmp_d_25.53, <retval>
.LVL16:
	.loc 1 23 11 view .LVU45
.LBE68:
	.loc 1 21 42 is_stmt 1 discriminator 3 view .LVU46
	.loc 1 21 30 discriminator 1 view .LVU47
	testb	$1, %cl	#, niters.41
	je	.L3	#,
	andl	$-2, %ecx	#, niters_vector_mult_vf.43
	addl	%ecx, %eax	# niters_vector_mult_vf.43,
.LVL17:
.L9:
.LBB69:
	.loc 1 22 9 view .LVU48
# src/kdtree.c:22:         FLOAT_TYPE dd = (p1[i] - p2[i]);
	.loc 1 22 20 is_stmt 0 view .LVU49
	vmovsd	(%rdx,%rax,8), %xmm1	# *_67, *_67
	vsubsd	(%rsi,%rax,8), %xmm1, %xmm1	# *_69, *_67, dd
.LVL18:
	.loc 1 23 9 is_stmt 1 view .LVU50
# src/kdtree.c:23:         d += dd*dd;
	.loc 1 23 16 is_stmt 0 view .LVU51
	vmulsd	%xmm1, %xmm1, %xmm1	# dd, dd, tmp148
.LVL19:
# src/kdtree.c:23:         d += dd*dd;
	.loc 1 23 11 view .LVU52
	vaddsd	%xmm1, %xmm0, %xmm0	# tmp148, <retval>, <retval>
.LVL20:
	.loc 1 23 11 view .LVU53
.LBE69:
	.loc 1 21 42 is_stmt 1 discriminator 3 view .LVU54
	.loc 1 21 30 discriminator 1 view .LVU55
	ret	
.LVL21:
	.p2align 4,,10
	.p2align 3
.L23:
	.loc 1 21 30 is_stmt 0 discriminator 1 view .LVU56
	vzeroupper
.LVL22:
.L3:
	.loc 1 21 30 discriminator 1 view .LVU57
.LBE71:
# src/kdtree.c:27: }
	.loc 1 27 1 view .LVU58
	ret	
.LVL23:
	.p2align 4,,10
	.p2align 3
.L11:
# src/kdtree.c:20:     FLOAT_TYPE d = 0;
	.loc 1 20 16 view .LVU59
	vxorpd	%xmm0, %xmm0, %xmm0	# <retval>
	.loc 1 25 2 is_stmt 1 view .LVU60
# src/kdtree.c:25: 	return d;
	.loc 1 25 9 is_stmt 0 view .LVU61
	ret	
.L12:
.LBB72:
# src/kdtree.c:21:     for(unsigned int i = 0; i<data_dims; ++i){
	.loc 1 21 22 view .LVU62
	xorl	%eax, %eax	#
.LBE72:
# src/kdtree.c:20:     FLOAT_TYPE d = 0;
	.loc 1 20 16 view .LVU63
	vxorpd	%xmm0, %xmm0, %xmm0	# <retval>
	jmp	.L5	#
	.cfi_endproc
.LFE54:
	.size	euclidean_distance, .-euclidean_distance
	.p2align 4
	.globl	swapHeapNode
	.type	swapHeapNode, @function
swapHeapNode:
.LVL24:
.LFB55:
	.loc 1 29 46 is_stmt 1 view -0
	.cfi_startproc
	.loc 1 29 46 is_stmt 0 view .LVU65
	endbr64	
	.loc 1 30 5 is_stmt 1 view .LVU66
	.loc 1 31 5 view .LVU67
.LVL25:
.LBB73:
.LBI73:
	.loc 2 26 1 view .LVU68
.LBB74:
	.loc 2 29 3 view .LVU69
# /usr/include/x86_64-linux-gnu/bits/string_fortified.h:29:   return __builtin___memcpy_chk (__dest, __src, __len,
	.loc 2 29 10 is_stmt 0 discriminator 1 view .LVU70
	vmovdqu	(%rdi), %xmm0	# MEM <uint128_t> [(char * {ref-all})a_2(D)], _3
.LVL26:
	.loc 2 29 10 discriminator 1 view .LVU71
.LBE74:
.LBE73:
	.loc 1 32 5 is_stmt 1 view .LVU72
.LBB75:
.LBI75:
	.loc 2 26 1 view .LVU73
.LBB76:
	.loc 2 29 3 view .LVU74
# /usr/include/x86_64-linux-gnu/bits/string_fortified.h:29:   return __builtin___memcpy_chk (__dest, __src, __len,
	.loc 2 29 10 is_stmt 0 discriminator 1 view .LVU75
	vmovdqu	(%rsi), %xmm1	# MEM <uint128_t> [(char * {ref-all})b_4(D)], tmp89
	vmovdqu	%xmm1, (%rdi)	# tmp89, MEM <uint128_t> [(char * {ref-all})a_2(D)]
.LVL27:
	.loc 2 29 10 discriminator 1 view .LVU76
.LBE76:
.LBE75:
	.loc 1 33 5 is_stmt 1 view .LVU77
.LBB77:
.LBI77:
	.loc 2 26 1 view .LVU78
.LBB78:
	.loc 2 29 3 view .LVU79
# /usr/include/x86_64-linux-gnu/bits/string_fortified.h:29:   return __builtin___memcpy_chk (__dest, __src, __len,
	.loc 2 29 10 is_stmt 0 discriminator 1 view .LVU80
	vmovdqu	%xmm0, (%rsi)	# _3, MEM <uint128_t> [(char * {ref-all})b_4(D)]
.LVL28:
	.loc 2 29 10 discriminator 1 view .LVU81
.LBE78:
.LBE77:
	.loc 1 34 5 is_stmt 1 view .LVU82
# src/kdtree.c:35: }
	.loc 1 35 1 is_stmt 0 view .LVU83
	ret	
	.cfi_endproc
.LFE55:
	.size	swapHeapNode, .-swapHeapNode
	.p2align 4
	.globl	swap_kd_node_ptrs
	.type	swap_kd_node_ptrs, @function
swap_kd_node_ptrs:
.LVL29:
.LFB56:
	.loc 1 37 50 is_stmt 1 view -0
	.cfi_startproc
	.loc 1 37 50 is_stmt 0 view .LVU85
	endbr64	
	.loc 1 38 5 is_stmt 1 view .LVU86
	.loc 1 39 5 view .LVU87
# src/kdtree.c:39:     tmp = *x;
	.loc 1 39 9 is_stmt 0 view .LVU88
	movq	(%rdi), %rax	# *x_3(D), tmp
.LVL30:
	.loc 1 40 5 is_stmt 1 view .LVU89
# src/kdtree.c:40:     *x = *y;
	.loc 1 40 10 is_stmt 0 view .LVU90
	movq	(%rsi), %rdx	# *y_5(D), _1
# src/kdtree.c:40:     *x = *y;
	.loc 1 40 8 view .LVU91
	movq	%rdx, (%rdi)	# _1, *x_3(D)
	.loc 1 41 5 is_stmt 1 view .LVU92
# src/kdtree.c:41:     *y = tmp;
	.loc 1 41 8 is_stmt 0 view .LVU93
	movq	%rax, (%rsi)	# tmp, *y_5(D)
# src/kdtree.c:45: }
	.loc 1 45 1 view .LVU94
	ret	
	.cfi_endproc
.LFE56:
	.size	swap_kd_node_ptrs, .-swap_kd_node_ptrs
	.p2align 4
	.globl	initializeKDnodes
	.type	initializeKDnodes, @function
initializeKDnodes:
.LVL31:
.LFB57:
	.loc 1 56 1 is_stmt 1 view -0
	.cfi_startproc
	.loc 1 56 1 is_stmt 0 view .LVU96
	endbr64	
	.loc 1 57 5 is_stmt 1 view .LVU97
.LBB79:
	.loc 1 57 9 view .LVU98
.LVL32:
	.loc 1 57 24 discriminator 1 view .LVU99
	testq	%rdx, %rdx	# n
	je	.L33	#,
# src/kdtree.c:59:         node_array[i].data = d + (i*data_dims);
	.loc 1 59 36 is_stmt 0 view .LVU100
	movl	data_dims(%rip), %ecx	# data_dims, data_dims
	movq	.LC1(%rip), %r8	#, tmp96
# src/kdtree.c:57:     for(idx_t i = 0; i < n; ++i)
	.loc 1 57 15 view .LVU101
	xorl	%eax, %eax	# i
# src/kdtree.c:63:         node_array[i].parent = NULL;
	.loc 1 63 30 view .LVU102
	vpxor	%xmm0, %xmm0, %xmm0	# tmp94
	salq	$3, %rcx	#, _20
.LVL33:
	.p2align 4,,10
	.p2align 3
.L28:
	.loc 1 59 9 is_stmt 1 view .LVU103
# src/kdtree.c:60:         node_array[i].array_idx = i;
	.loc 1 60 33 is_stmt 0 view .LVU104
	movq	%rax, 16(%rdi)	# i, MEM[(long unsigned int *)_32 + 16B]
# src/kdtree.c:57:     for(idx_t i = 0; i < n; ++i)
	.loc 1 57 29 discriminator 3 view .LVU105
	addq	$1, %rax	#, i
.LVL34:
# src/kdtree.c:57:     for(idx_t i = 0; i < n; ++i)
	.loc 1 57 24 discriminator 1 view .LVU106
	addq	$48, %rdi	#, ivtmp.84
.LVL35:
# src/kdtree.c:59:         node_array[i].data = d + (i*data_dims);
	.loc 1 59 28 view .LVU107
	movq	%rsi, -40(%rdi)	# ivtmp.83, MEM[(double * *)_32 + 8B]
	.loc 1 60 9 is_stmt 1 view .LVU108
	.loc 1 61 9 view .LVU109
	.loc 1 62 9 view .LVU110
# src/kdtree.c:57:     for(idx_t i = 0; i < n; ++i)
	.loc 1 57 24 is_stmt 0 discriminator 1 view .LVU111
	addq	%rcx, %rsi	# _20, ivtmp.83
# src/kdtree.c:62:         node_array[i].rch = NULL;
	.loc 1 62 27 view .LVU112
	movq	$0, -8(%rdi)	#, MEM[(struct kd_node * *)_32 + 40B]
	.loc 1 63 9 is_stmt 1 view .LVU113
# src/kdtree.c:63:         node_array[i].parent = NULL;
	.loc 1 63 30 is_stmt 0 view .LVU114
	vmovdqu	%xmm0, -24(%rdi)	# tmp94, MEM <vector(2) long unsigned int> [(struct kd_node * *)_32 + 24B]
	.loc 1 64 9 is_stmt 1 view .LVU115
	.loc 1 65 9 view .LVU116
# src/kdtree.c:64:         node_array[i].level = -1;
	.loc 1 64 29 is_stmt 0 view .LVU117
	movq	%r8, -48(%rdi)	# tmp96, MEM <vector(2) int> [(int *)_32]
	.loc 1 57 29 is_stmt 1 discriminator 3 view .LVU118
.LVL36:
	.loc 1 57 24 discriminator 1 view .LVU119
	cmpq	%rax, %rdx	# i, n
	jne	.L28	#,
.LVL37:
.L33:
	.loc 1 57 24 is_stmt 0 discriminator 1 view .LVU120
.LBE79:
# src/kdtree.c:67: }
	.loc 1 67 1 view .LVU121
	ret	
	.cfi_endproc
.LFE57:
	.size	initializeKDnodes, .-initializeKDnodes
	.p2align 4
	.globl	initializePTRS
	.type	initializePTRS, @function
initializePTRS:
.LVL38:
.LFB58:
	.loc 1 70 1 is_stmt 1 view -0
	.cfi_startproc
	.loc 1 70 1 is_stmt 0 view .LVU123
	endbr64	
	.loc 1 71 5 is_stmt 1 view .LVU124
.LBB80:
	.loc 1 71 9 view .LVU125
.LVL39:
	.loc 1 71 24 discriminator 1 view .LVU126
.LBE80:
# src/kdtree.c:70: {
	.loc 1 70 1 is_stmt 0 view .LVU127
	movq	%rdx, %rcx	# tmp140, n
.LBB81:
# src/kdtree.c:71:     for(idx_t i = 0; i < n; ++i)
	.loc 1 71 24 discriminator 1 view .LVU128
	testq	%rdx, %rdx	# n
	je	.L48	#,
	leaq	-1(%rdx), %rax	#, tmp113
	cmpq	$2, %rax	#, tmp113
	jbe	.L39	#,
	shrq	$2, %rdx	#, bnd.90
.LVL40:
	.loc 1 71 24 discriminator 1 view .LVU129
	vmovq	%rsi, %xmm5	# node_array, node_array
	movq	%rdi, %rax	# node_ptr_array, ivtmp.102
	vmovdqa	.LC3(%rip), %ymm1	#, vect_vec_iv_.93
	salq	$5, %rdx	#, tmp116
	vpbroadcastq	%xmm5, %ymm4	# node_array, vect_cst__36
	vpbroadcastq	.LC5(%rip), %ymm3	#, tmp118
	addq	%rdi, %rdx	# node_ptr_array, _50
.LVL41:
	.p2align 4,,10
	.p2align 3
.L37:
	.loc 1 71 24 discriminator 1 view .LVU130
	vmovdqa	%ymm1, %ymm2	# vect_vec_iv_.93, vect_vec_iv_.93
	addq	$32, %rax	#, ivtmp.102
	vpaddq	%ymm3, %ymm1, %ymm1	# tmp118, vect_vec_iv_.93, vect_vec_iv_.93
	.loc 1 73 9 is_stmt 1 view .LVU131
# src/kdtree.c:73:         node_ptr_array[i] = node_array + i;
	.loc 1 73 40 is_stmt 0 view .LVU132
	vpsllq	$1, %ymm2, %ymm0	#, vect_vec_iv_.93, tmp121
	vpaddq	%ymm2, %ymm0, %ymm0	# vect_vec_iv_.93, tmp121, vect__1.94
	vpsllq	$4, %ymm0, %ymm0	#, vect__1.94, tmp123
	vpaddq	%ymm4, %ymm0, %ymm0	# vect_cst__36, tmp123, vect__4.95
# src/kdtree.c:73:         node_ptr_array[i] = node_array + i;
	.loc 1 73 27 view .LVU133
	vmovdqu	%ymm0, -32(%rax)	# vect__4.95, MEM <vector(4) long unsigned int> [(struct kd_node * *)_18]
	.loc 1 71 29 is_stmt 1 discriminator 3 view .LVU134
	.loc 1 71 24 discriminator 1 view .LVU135
	cmpq	%rdx, %rax	# _50, ivtmp.102
	jne	.L37	#,
	testb	$3, %cl	#, n
	je	.L47	#,
	movq	%rcx, %rax	# n, niters_vector_mult_vf.91
	andq	$-4, %rax	#, niters_vector_mult_vf.91
	vzeroupper
.LVL42:
.L36:
	.loc 1 73 9 view .LVU136
# src/kdtree.c:73:         node_ptr_array[i] = node_array + i;
	.loc 1 73 40 is_stmt 0 view .LVU137
	leaq	(%rax,%rax,2), %rdx	#, tmp128
# src/kdtree.c:73:         node_ptr_array[i] = node_array + i;
	.loc 1 73 23 view .LVU138
	leaq	0(,%rax,8), %r8	#, _3
# src/kdtree.c:73:         node_ptr_array[i] = node_array + i;
	.loc 1 73 40 view .LVU139
	salq	$4, %rdx	#, tmp129
	leaq	(%rsi,%rdx), %r9	#, tmp130
	movq	%r9, (%rdi,%rax,8)	# tmp130, *_4
	.loc 1 71 29 is_stmt 1 discriminator 3 view .LVU140
	.loc 1 71 24 discriminator 1 view .LVU141
# src/kdtree.c:71:     for(idx_t i = 0; i < n; ++i)
	.loc 1 71 29 is_stmt 0 discriminator 3 view .LVU142
	leaq	1(%rax), %r9	#, i
# src/kdtree.c:71:     for(idx_t i = 0; i < n; ++i)
	.loc 1 71 24 discriminator 1 view .LVU143
	cmpq	%rcx, %r9	# n, i
	jnb	.L48	#,
	.loc 1 73 9 is_stmt 1 view .LVU144
# src/kdtree.c:73:         node_ptr_array[i] = node_array + i;
	.loc 1 73 40 is_stmt 0 view .LVU145
	leaq	48(%rsi,%rdx), %r9	#, tmp133
# src/kdtree.c:71:     for(idx_t i = 0; i < n; ++i)
	.loc 1 71 29 discriminator 3 view .LVU146
	addq	$2, %rax	#, i
# src/kdtree.c:73:         node_ptr_array[i] = node_array + i;
	.loc 1 73 40 view .LVU147
	movq	%r9, 8(%rdi,%r8)	# tmp133, *_45
	.loc 1 71 29 is_stmt 1 discriminator 3 view .LVU148
	.loc 1 71 24 discriminator 1 view .LVU149
	cmpq	%rcx, %rax	# n, i
	jnb	.L48	#,
	.loc 1 73 9 view .LVU150
# src/kdtree.c:73:         node_ptr_array[i] = node_array + i;
	.loc 1 73 40 is_stmt 0 view .LVU151
	leaq	96(%rsi,%rdx), %rax	#, tmp136
	movq	%rax, 16(%rdi,%r8)	# tmp136, *_5
	.loc 1 71 29 is_stmt 1 discriminator 3 view .LVU152
	.loc 1 71 24 discriminator 1 view .LVU153
.LBE81:
# src/kdtree.c:75: }
	.loc 1 75 1 is_stmt 0 view .LVU154
	ret	
.LVL43:
	.p2align 4,,10
	.p2align 3
.L47:
	.loc 1 75 1 view .LVU155
	vzeroupper
.LVL44:
.L48:
	.loc 1 75 1 view .LVU156
	ret	
.LVL45:
.L39:
.LBB82:
# src/kdtree.c:71:     for(idx_t i = 0; i < n; ++i)
	.loc 1 71 15 view .LVU157
	xorl	%eax, %eax	# niters_vector_mult_vf.91
	jmp	.L36	#
.LBE82:
	.cfi_endproc
.LFE58:
	.size	initializePTRS, .-initializePTRS
	.p2align 4
	.globl	cmpKDnodes
	.type	cmpKDnodes, @function
cmpKDnodes:
.LVL46:
.LFB59:
	.loc 1 77 48 is_stmt 1 view -0
	.cfi_startproc
	.loc 1 77 48 is_stmt 0 view .LVU159
	endbr64	
	.loc 1 79 5 is_stmt 1 view .LVU160
# src/kdtree.c:79:     FLOAT_TYPE res = a->data[var] - b->data[var];
	.loc 1 79 44 is_stmt 0 view .LVU161
	movq	8(%rsi), %rax	# b_13(D)->data, b_13(D)->data
# src/kdtree.c:79:     FLOAT_TYPE res = a->data[var] - b->data[var];
	.loc 1 79 29 view .LVU162
	movq	8(%rdi), %rcx	# a_11(D)->data, a_11(D)->data
	movslq	%edx, %rdx	# tmp107, var
.LVL47:
	.loc 1 80 5 is_stmt 1 view .LVU163
# src/kdtree.c:80:     return (res > 0);
	.loc 1 80 17 is_stmt 0 view .LVU164
	vxorpd	%xmm1, %xmm1, %xmm1	# tmp104
# src/kdtree.c:79:     FLOAT_TYPE res = a->data[var] - b->data[var];
	.loc 1 79 16 view .LVU165
	vmovsd	(%rcx,%rdx,8), %xmm0	# *_4, *_4
	vsubsd	(%rax,%rdx,8), %xmm0, %xmm0	# *_7, *_4, res
# src/kdtree.c:80:     return (res > 0);
	.loc 1 80 17 view .LVU166
	xorl	%eax, %eax	# tmp102
	vcomisd	%xmm1, %xmm0	# tmp104, res
	seta	%al	#, tmp102
# src/kdtree.c:81: }
	.loc 1 81 1 view .LVU167
	ret	
	.cfi_endproc
.LFE59:
	.size	cmpKDnodes, .-cmpKDnodes
	.section	.rodata.str1.1,"aMS",@progbits,1
.LC6:
	.string	"Node %p:\n"
.LC7:
	.string	"\t array_idx: %lu\n"
.LC8:
	.string	"\t data: "
.LC9:
	.string	" %f "
.LC10:
	.string	"\t parent: %p\n"
.LC11:
	.string	"\t lch: %p\n"
.LC12:
	.string	"\t rch: %p\n"
.LC13:
	.string	"\t level: %d\n"
.LC14:
	.string	"\t split_var: %d\n"
	.text
	.p2align 4
	.globl	printKDnode
	.type	printKDnode, @function
printKDnode:
.LVL48:
.LFB60:
	.loc 1 84 1 is_stmt 1 view -0
	.cfi_startproc
	.loc 1 84 1 is_stmt 0 view .LVU169
	endbr64	
	.loc 1 85 5 is_stmt 1 view .LVU170
.LVL49:
.LBB83:
.LBI83:
	.file 3 "/usr/include/x86_64-linux-gnu/bits/stdio2.h"
	.loc 3 84 1 view .LVU171
.LBB84:
	.loc 3 86 3 view .LVU172
.LBE84:
.LBE83:
# src/kdtree.c:84: {
	.loc 1 84 1 is_stmt 0 view .LVU173
	pushq	%r12	#
	.cfi_def_cfa_offset 16
	.cfi_offset 12, -16
.LBB88:
.LBB85:
# /usr/include/x86_64-linux-gnu/bits/stdio2.h:86:   return __printf_chk (__USE_FORTIFY_LEVEL - 1, __fmt, __va_arg_pack ());
	.loc 3 86 10 view .LVU174
	movq	%rdi, %rdx	# node,
	leaq	.LC6(%rip), %rsi	#, tmp97
	xorl	%eax, %eax	#
.LBE85:
.LBE88:
# src/kdtree.c:84: {
	.loc 1 84 1 view .LVU175
	pushq	%rbp	#
	.cfi_def_cfa_offset 24
	.cfi_offset 6, -24
	movq	%rdi, %rbp	# tmp115, node
.LBB89:
.LBB86:
# /usr/include/x86_64-linux-gnu/bits/stdio2.h:86:   return __printf_chk (__USE_FORTIFY_LEVEL - 1, __fmt, __va_arg_pack ());
	.loc 3 86 10 view .LVU176
	movl	$1, %edi	#,
.LVL50:
	.loc 3 86 10 view .LVU177
.LBE86:
.LBE89:
# src/kdtree.c:84: {
	.loc 1 84 1 view .LVU178
	pushq	%rbx	#
	.cfi_def_cfa_offset 32
	.cfi_offset 3, -32
.LBB90:
.LBB87:
# /usr/include/x86_64-linux-gnu/bits/stdio2.h:86:   return __printf_chk (__USE_FORTIFY_LEVEL - 1, __fmt, __va_arg_pack ());
	.loc 3 86 10 view .LVU179
	call	__printf_chk@PLT	#
.LVL51:
	.loc 3 86 10 view .LVU180
.LBE87:
.LBE90:
	.loc 1 86 5 is_stmt 1 view .LVU181
.LBB91:
.LBI91:
	.loc 3 84 1 view .LVU182
.LBB92:
	.loc 3 86 3 view .LVU183
# /usr/include/x86_64-linux-gnu/bits/stdio2.h:86:   return __printf_chk (__USE_FORTIFY_LEVEL - 1, __fmt, __va_arg_pack ());
	.loc 3 86 10 is_stmt 0 view .LVU184
	movq	16(%rbp), %rdx	# node_16(D)->array_idx, node_16(D)->array_idx
	movl	$1, %edi	#,
	xorl	%eax, %eax	#
	leaq	.LC7(%rip), %rsi	#, tmp99
	call	__printf_chk@PLT	#
.LVL52:
	.loc 3 86 10 view .LVU185
.LBE92:
.LBE91:
	.loc 1 87 5 is_stmt 1 view .LVU186
.LBB93:
.LBI93:
	.loc 3 84 1 view .LVU187
.LBB94:
	.loc 3 86 3 view .LVU188
# /usr/include/x86_64-linux-gnu/bits/stdio2.h:86:   return __printf_chk (__USE_FORTIFY_LEVEL - 1, __fmt, __va_arg_pack ());
	.loc 3 86 10 is_stmt 0 view .LVU189
	leaq	.LC8(%rip), %rsi	#, tmp100
	movl	$1, %edi	#,
	xorl	%eax, %eax	#
	call	__printf_chk@PLT	#
.LVL53:
	.loc 3 86 10 view .LVU190
.LBE94:
.LBE93:
	.loc 1 88 5 is_stmt 1 view .LVU191
.LBB95:
	.loc 1 88 9 view .LVU192
	.loc 1 88 28 discriminator 2 view .LVU193
	movl	data_dims(%rip), %eax	# data_dims,
	testl	%eax, %eax	#
	je	.L51	#,
	xorl	%ebx, %ebx	# ivtmp.114
	leaq	.LC9(%rip), %r12	#, tmp114
.LVL54:
	.p2align 4,,10
	.p2align 3
.L52:
	.loc 1 88 45 discriminator 3 view .LVU194
.LBB96:
.LBI96:
	.loc 3 84 1 view .LVU195
.LBB97:
	.loc 3 86 3 view .LVU196
.LBE97:
.LBE96:
# src/kdtree.c:88:     for(unsigned int i=0; i<data_dims; ++i) printf(" %f ",node->data[i]);
	.loc 1 88 45 is_stmt 0 discriminator 3 view .LVU197
	movq	8(%rbp), %rax	# node_16(D)->data, node_16(D)->data
.LBB100:
.LBB98:
# /usr/include/x86_64-linux-gnu/bits/stdio2.h:86:   return __printf_chk (__USE_FORTIFY_LEVEL - 1, __fmt, __va_arg_pack ());
	.loc 3 86 10 view .LVU198
	movq	%r12, %rsi	# tmp114,
	movl	$1, %edi	#,
	vmovsd	(%rax,%rbx,8), %xmm0	# *_5, *_5
	movl	$1, %eax	#,
.LBE98:
.LBE100:
# src/kdtree.c:88:     for(unsigned int i=0; i<data_dims; ++i) printf(" %f ",node->data[i]);
	.loc 1 88 28 discriminator 2 view .LVU199
	addq	$1, %rbx	#, ivtmp.114
.LVL55:
.LBB101:
.LBB99:
# /usr/include/x86_64-linux-gnu/bits/stdio2.h:86:   return __printf_chk (__USE_FORTIFY_LEVEL - 1, __fmt, __va_arg_pack ());
	.loc 3 86 10 view .LVU200
	call	__printf_chk@PLT	#
.LVL56:
	.loc 3 86 10 view .LVU201
.LBE99:
.LBE101:
	.loc 1 88 40 is_stmt 1 discriminator 3 view .LVU202
	.loc 1 88 28 discriminator 2 view .LVU203
	cmpl	data_dims(%rip), %ebx	# data_dims, ivtmp.114
	jb	.L52	#,
.L51:
	.loc 1 88 28 is_stmt 0 discriminator 2 view .LVU204
.LBE95:
	.loc 1 89 5 is_stmt 1 view .LVU205
.LVL57:
.LBB102:
.LBI102:
	.loc 3 84 1 view .LVU206
.LBB103:
	.loc 3 86 3 view .LVU207
# /usr/include/x86_64-linux-gnu/bits/stdio2.h:86:   return __printf_chk (__USE_FORTIFY_LEVEL - 1, __fmt, __va_arg_pack ());
	.loc 3 86 10 is_stmt 0 view .LVU208
	movl	$10, %edi	#,
	call	putchar@PLT	#
.LVL58:
	.loc 3 86 10 view .LVU209
.LBE103:
.LBE102:
	.loc 1 90 5 is_stmt 1 view .LVU210
.LBB104:
.LBI104:
	.loc 3 84 1 view .LVU211
.LBB105:
	.loc 3 86 3 view .LVU212
# /usr/include/x86_64-linux-gnu/bits/stdio2.h:86:   return __printf_chk (__USE_FORTIFY_LEVEL - 1, __fmt, __va_arg_pack ());
	.loc 3 86 10 is_stmt 0 view .LVU213
	movq	24(%rbp), %rdx	# node_16(D)->parent, node_16(D)->parent
	movl	$1, %edi	#,
	xorl	%eax, %eax	#
	leaq	.LC10(%rip), %rsi	#, tmp105
	call	__printf_chk@PLT	#
.LVL59:
	.loc 3 86 10 view .LVU214
.LBE105:
.LBE104:
	.loc 1 91 5 is_stmt 1 view .LVU215
.LBB106:
.LBI106:
	.loc 3 84 1 view .LVU216
.LBB107:
	.loc 3 86 3 view .LVU217
# /usr/include/x86_64-linux-gnu/bits/stdio2.h:86:   return __printf_chk (__USE_FORTIFY_LEVEL - 1, __fmt, __va_arg_pack ());
	.loc 3 86 10 is_stmt 0 view .LVU218
	movq	32(%rbp), %rdx	# node_16(D)->lch, node_16(D)->lch
	movl	$1, %edi	#,
	xorl	%eax, %eax	#
	leaq	.LC11(%rip), %rsi	#, tmp107
	call	__printf_chk@PLT	#
.LVL60:
	.loc 3 86 10 view .LVU219
.LBE107:
.LBE106:
	.loc 1 92 5 is_stmt 1 view .LVU220
.LBB108:
.LBI108:
	.loc 3 84 1 view .LVU221
.LBB109:
	.loc 3 86 3 view .LVU222
# /usr/include/x86_64-linux-gnu/bits/stdio2.h:86:   return __printf_chk (__USE_FORTIFY_LEVEL - 1, __fmt, __va_arg_pack ());
	.loc 3 86 10 is_stmt 0 view .LVU223
	movq	40(%rbp), %rdx	# node_16(D)->rch, node_16(D)->rch
	movl	$1, %edi	#,
	xorl	%eax, %eax	#
	leaq	.LC12(%rip), %rsi	#, tmp109
	call	__printf_chk@PLT	#
.LVL61:
	.loc 3 86 10 view .LVU224
.LBE109:
.LBE108:
	.loc 1 93 5 is_stmt 1 view .LVU225
.LBB110:
.LBI110:
	.loc 3 84 1 view .LVU226
.LBB111:
	.loc 3 86 3 view .LVU227
# /usr/include/x86_64-linux-gnu/bits/stdio2.h:86:   return __printf_chk (__USE_FORTIFY_LEVEL - 1, __fmt, __va_arg_pack ());
	.loc 3 86 10 is_stmt 0 view .LVU228
	movl	0(%rbp), %edx	# node_16(D)->level, node_16(D)->level
	movl	$1, %edi	#,
	xorl	%eax, %eax	#
	leaq	.LC13(%rip), %rsi	#, tmp111
	call	__printf_chk@PLT	#
.LVL62:
	.loc 3 86 10 view .LVU229
.LBE111:
.LBE110:
	.loc 1 94 5 is_stmt 1 view .LVU230
.LBB112:
.LBI112:
	.loc 3 84 1 view .LVU231
.LBB113:
	.loc 3 86 3 view .LVU232
# /usr/include/x86_64-linux-gnu/bits/stdio2.h:86:   return __printf_chk (__USE_FORTIFY_LEVEL - 1, __fmt, __va_arg_pack ());
	.loc 3 86 10 is_stmt 0 view .LVU233
	movl	4(%rbp), %edx	# node_16(D)->split_var, node_16(D)->split_var
	movl	$1, %edi	#,
	xorl	%eax, %eax	#
	leaq	.LC14(%rip), %rsi	#, tmp113
	call	__printf_chk@PLT	#
.LVL63:
	.loc 3 86 10 view .LVU234
.LBE113:
.LBE112:
	.loc 1 95 5 is_stmt 1 view .LVU235
.LBB114:
.LBI114:
	.loc 3 84 1 view .LVU236
.LBB115:
	.loc 3 86 3 view .LVU237
.LBE115:
.LBE114:
# src/kdtree.c:96: }
	.loc 1 96 1 is_stmt 0 view .LVU238
	popq	%rbx	#
	.cfi_def_cfa_offset 24
.LBB118:
.LBB116:
# /usr/include/x86_64-linux-gnu/bits/stdio2.h:86:   return __printf_chk (__USE_FORTIFY_LEVEL - 1, __fmt, __va_arg_pack ());
	.loc 3 86 10 view .LVU239
	movl	$10, %edi	#,
.LBE116:
.LBE118:
# src/kdtree.c:96: }
	.loc 1 96 1 view .LVU240
	popq	%rbp	#
	.cfi_def_cfa_offset 16
.LVL64:
	.loc 1 96 1 view .LVU241
	popq	%r12	#
	.cfi_def_cfa_offset 8
.LBB119:
.LBB117:
# /usr/include/x86_64-linux-gnu/bits/stdio2.h:86:   return __printf_chk (__USE_FORTIFY_LEVEL - 1, __fmt, __va_arg_pack ());
	.loc 3 86 10 view .LVU242
	jmp	putchar@PLT	#
.LVL65:
.LBE117:
.LBE119:
	.cfi_endproc
.LFE60:
	.size	printKDnode, .-printKDnode
	.p2align 4
	.globl	partition
	.type	partition, @function
partition:
.LVL66:
.LFB61:
	.loc 1 101 1 is_stmt 1 view -0
	.cfi_startproc
	.loc 1 101 1 is_stmt 0 view .LVU244
	endbr64	
# src/kdtree.c:101: {
	.loc 1 101 1 view .LVU245
	movl	%esi, %eax	# tmp144, low
	.loc 1 102 5 is_stmt 1 view .LVU246
# src/kdtree.c:102:     kd_node* pivot = arr[high];
	.loc 1 102 25 is_stmt 0 view .LVU247
	movslq	%edx, %rsi	# high, high
.LVL67:
	.loc 1 102 25 view .LVU248
	leaq	(%rdi,%rsi,8), %r10	#, _3
# src/kdtree.c:104:     int i = (low - 1);
	.loc 1 104 9 view .LVU249
	leal	-1(%rax), %esi	#, i
# src/kdtree.c:102:     kd_node* pivot = arr[high];
	.loc 1 102 14 view .LVU250
	movq	(%r10), %r8	# *_3, pivot
.LVL68:
	.loc 1 104 5 is_stmt 1 view .LVU251
	.loc 1 105 5 view .LVU252
.LBB127:
	.loc 1 105 10 view .LVU253
	.loc 1 105 25 discriminator 1 view .LVU254
	cmpl	%eax, %edx	# low, high
	jle	.L56	#,
.LBB128:
.LBB129:
# src/kdtree.c:79:     FLOAT_TYPE res = a->data[var] - b->data[var];
	.loc 1 79 44 is_stmt 0 view .LVU255
	movq	8(%r8), %r8	# pivot_24->data, pivot_24->data
.LVL69:
# src/kdtree.c:79:     FLOAT_TYPE res = a->data[var] - b->data[var];
	.loc 1 79 29 view .LVU256
	movslq	%ecx, %rcx	# split_var, split_var
	.loc 1 79 29 view .LVU257
	subl	%eax, %edx	# low, tmp130
.LVL70:
	.loc 1 79 29 view .LVU258
.LBE129:
.LBE128:
# src/kdtree.c:106:         if (!cmpKDnodes(arr[j],pivot,split_var)) {
	.loc 1 106 12 discriminator 1 view .LVU259
	vxorpd	%xmm1, %xmm1, %xmm1	# tmp136
.LBB132:
.LBB130:
# src/kdtree.c:79:     FLOAT_TYPE res = a->data[var] - b->data[var];
	.loc 1 79 29 view .LVU260
	leaq	0(,%rcx,8), %r9	#, _33
# src/kdtree.c:79:     FLOAT_TYPE res = a->data[var] - b->data[var];
	.loc 1 79 44 view .LVU261
	vmovsd	(%r8,%rcx,8), %xmm2	# *_37, _38
	movslq	%eax, %r8	# low, _51
	addq	%r8, %rdx	# _51, tmp131
	leaq	(%rdi,%r8,8), %rcx	#, ivtmp.126
.LVL71:
	.loc 1 79 44 view .LVU262
	leaq	(%rdi,%rdx,8), %r8	#, _91
.LVL72:
	.p2align 4,,10
	.p2align 3
.L58:
	.loc 1 79 44 view .LVU263
.LBE130:
.LBE132:
	.loc 1 106 9 is_stmt 1 view .LVU264
# src/kdtree.c:106:         if (!cmpKDnodes(arr[j],pivot,split_var)) {
	.loc 1 106 14 is_stmt 0 view .LVU265
	movq	(%rcx), %rax	# MEM[(struct kd_node * *)_5], _75
.LVL73:
.LBB133:
.LBI128:
	.loc 1 77 5 is_stmt 1 view .LVU266
.LBB131:
	.loc 1 79 5 view .LVU267
	.loc 1 80 5 view .LVU268
# src/kdtree.c:79:     FLOAT_TYPE res = a->data[var] - b->data[var];
	.loc 1 79 29 is_stmt 0 view .LVU269
	movq	8(%rax), %rdx	# _75->data, _75->data
# src/kdtree.c:79:     FLOAT_TYPE res = a->data[var] - b->data[var];
	.loc 1 79 16 view .LVU270
	vmovsd	(%rdx,%r9), %xmm0	# *_73, *_73
	vsubsd	%xmm2, %xmm0, %xmm0	# _38, *_73, res
.LBE131:
.LBE133:
# src/kdtree.c:106:         if (!cmpKDnodes(arr[j],pivot,split_var)) {
	.loc 1 106 12 discriminator 1 view .LVU271
	vcomisd	%xmm1, %xmm0	# tmp136, res
	ja	.L57	#,
# src/kdtree.c:107:             i++;
	.loc 1 107 14 view .LVU272
	addl	$1, %esi	#, i
.LVL74:
	.loc 1 107 13 is_stmt 1 view .LVU273
	.loc 1 108 13 view .LVU274
# src/kdtree.c:108:             swap_kd_node_ptrs(arr + i, arr + j);
	.loc 1 108 35 is_stmt 0 view .LVU275
	movslq	%esi, %rdx	# i, i
# src/kdtree.c:108:             swap_kd_node_ptrs(arr + i, arr + j);
	.loc 1 108 13 view .LVU276
	leaq	(%rdi,%rdx,8), %rdx	#, _67
.LVL75:
.LBB134:
.LBI134:
	.loc 1 37 6 is_stmt 1 view .LVU277
.LBB135:
	.loc 1 38 5 view .LVU278
	.loc 1 39 5 view .LVU279
# src/kdtree.c:39:     tmp = *x;
	.loc 1 39 9 is_stmt 0 view .LVU280
	movq	(%rdx), %r11	# *_67, tmp
.LVL76:
	.loc 1 40 5 is_stmt 1 view .LVU281
# src/kdtree.c:40:     *x = *y;
	.loc 1 40 8 is_stmt 0 view .LVU282
	movq	%rax, (%rdx)	# _75, *_67
	.loc 1 41 5 is_stmt 1 view .LVU283
# src/kdtree.c:41:     *y = tmp;
	.loc 1 41 8 is_stmt 0 view .LVU284
	movq	%r11, (%rcx)	# tmp, MEM[(struct kd_node * *)_5]
.LVL77:
.L57:
	.loc 1 41 8 view .LVU285
.LBE135:
.LBE134:
	.loc 1 105 39 is_stmt 1 discriminator 2 view .LVU286
	.loc 1 105 25 discriminator 1 view .LVU287
	addq	$8, %rcx	#, ivtmp.126
	cmpq	%r8, %rcx	# _91, ivtmp.126
	jne	.L58	#,
.LBE127:
.LBB136:
.LBB137:
# src/kdtree.c:40:     *x = *y;
	.loc 1 40 10 is_stmt 0 view .LVU288
	movq	(%r10), %r8	# *_3, pivot
.LBE137:
.LBE136:
# src/kdtree.c:112:     return (i + 1);
	.loc 1 112 15 view .LVU289
	leal	1(%rsi), %eax	#, <retval>
.LVL78:
.L56:
	.loc 1 111 5 is_stmt 1 view .LVU290
# src/kdtree.c:111:     swap_kd_node_ptrs(arr + i + 1, arr + high);
	.loc 1 111 31 is_stmt 0 view .LVU291
	movslq	%esi, %rsi	# i, i
# src/kdtree.c:111:     swap_kd_node_ptrs(arr + i + 1, arr + high);
	.loc 1 111 5 view .LVU292
	leaq	8(%rdi,%rsi,8), %rdx	#, _15
.LVL79:
.LBB139:
.LBI136:
	.loc 1 37 6 is_stmt 1 view .LVU293
.LBB138:
	.loc 1 38 5 view .LVU294
	.loc 1 39 5 view .LVU295
# src/kdtree.c:39:     tmp = *x;
	.loc 1 39 9 is_stmt 0 view .LVU296
	movq	(%rdx), %rcx	# *_15, tmp
.LVL80:
	.loc 1 40 5 is_stmt 1 view .LVU297
# src/kdtree.c:40:     *x = *y;
	.loc 1 40 8 is_stmt 0 view .LVU298
	movq	%r8, (%rdx)	# pivot, *_15
	.loc 1 41 5 is_stmt 1 view .LVU299
# src/kdtree.c:41:     *y = tmp;
	.loc 1 41 8 is_stmt 0 view .LVU300
	movq	%rcx, (%r10)	# tmp, *_3
.LVL81:
	.loc 1 41 8 view .LVU301
.LBE138:
.LBE139:
	.loc 1 112 5 is_stmt 1 view .LVU302
# src/kdtree.c:113: }
	.loc 1 113 1 is_stmt 0 view .LVU303
	ret	
	.cfi_endproc
.LFE61:
	.size	partition, .-partition
	.p2align 4
	.globl	medianOfNodes
	.type	medianOfNodes, @function
medianOfNodes:
.LVL82:
.LFB62:
	.loc 1 117 1 is_stmt 1 view -0
	.cfi_startproc
	.loc 1 117 1 is_stmt 0 view .LVU305
	endbr64	
	.loc 1 119 5 is_stmt 1 view .LVU306
.LVL83:
	.loc 1 129 5 view .LVU307
# src/kdtree.c:129:     if(left == right) return left;
	.loc 1 129 30 is_stmt 0 discriminator 1 view .LVU308
	movl	%edx, %r8d	# right, <retval>
# src/kdtree.c:129:     if(left == right) return left;
	.loc 1 129 7 view .LVU309
	cmpl	%esi, %edx	# left, right
	je	.L80	#,
# src/kdtree.c:130:     if(left == (right - 1)){
	.loc 1 130 23 view .LVU310
	leal	-1(%rdx), %eax	#, tmp146
	movq	%rdi, %r9	# tmp182, a
	movl	%esi, %r10d	# tmp183, left
	movl	%edx, %r11d	# tmp184, right
	.loc 1 130 5 is_stmt 1 view .LVU311
# src/kdtree.c:130:     if(left == (right - 1)){
	.loc 1 130 7 is_stmt 0 view .LVU312
	cmpl	%esi, %eax	# left, tmp146
	je	.L63	#,
	.loc 1 139 17 is_stmt 1 view .LVU313
	cmpl	%esi, %edx	# left, right
	jl	.L83	#,
# src/kdtree.c:117: {
	.loc 1 117 1 is_stmt 0 view .LVU314
	pushq	%r14	#
	.cfi_def_cfa_offset 16
	.cfi_offset 14, -16
.LBB157:
.LBB158:
.LBB159:
.LBB160:
.LBB161:
.LBB162:
# src/kdtree.c:79:     FLOAT_TYPE res = a->data[var] - b->data[var];
	.loc 1 79 29 view .LVU315
	movslq	%ecx, %rdi	# split_var, split_var
.LVL84:
	.loc 1 79 29 view .LVU316
.LBE162:
.LBE161:
# src/kdtree.c:106:         if (!cmpKDnodes(arr[j],pivot,split_var)) {
	.loc 1 106 12 discriminator 1 view .LVU317
	vxorpd	%xmm2, %xmm2, %xmm2	# tmp181
.LBE160:
.LBE159:
.LBE158:
.LBE157:
# src/kdtree.c:117: {
	.loc 1 117 1 view .LVU318
	pushq	%r13	#
	.cfi_def_cfa_offset 24
	.cfi_offset 13, -24
.LBB184:
.LBB181:
.LBB178:
.LBB171:
.LBB166:
.LBB163:
# src/kdtree.c:79:     FLOAT_TYPE res = a->data[var] - b->data[var];
	.loc 1 79 29 view .LVU319
	salq	$3, %rdi	#, _55
	movslq	%edx, %r13	# right, right
.LBE163:
.LBE166:
.LBE171:
.LBE178:
.LBE181:
.LBE184:
# src/kdtree.c:117: {
	.loc 1 117 1 view .LVU320
	pushq	%r12	#
	.cfi_def_cfa_offset 32
	.cfi_offset 12, -32
	leal	-1(%rsi), %r12d	#, i
	pushq	%rbp	#
	.cfi_def_cfa_offset 40
	.cfi_offset 6, -40
# src/kdtree.c:119:     int k = left + ((right - left + 1)/2); 
	.loc 1 119 28 view .LVU321
	movl	%edx, %ebp	# right, tmp157
	subl	%esi, %ebp	# left, tmp157
# src/kdtree.c:117: {
	.loc 1 117 1 view .LVU322
	pushq	%rbx	#
	.cfi_def_cfa_offset 48
	.cfi_offset 3, -48
# src/kdtree.c:119:     int k = left + ((right - left + 1)/2); 
	.loc 1 119 35 view .LVU323
	addl	$1, %ebp	#, tmp158
# src/kdtree.c:119:     int k = left + ((right - left + 1)/2); 
	.loc 1 119 39 view .LVU324
	sarl	%ebp	# tmp159
# src/kdtree.c:119:     int k = left + ((right - left + 1)/2); 
	.loc 1 119 9 view .LVU325
	addl	%esi, %ebp	# left, k
.LVL85:
	.p2align 4,,10
	.p2align 3
.L73:
.LBB185:
	.loc 1 143 9 is_stmt 1 view .LVU326
.LBB182:
.LBI158:
	.loc 1 100 5 view .LVU327
.LBB179:
	.loc 1 102 5 view .LVU328
# src/kdtree.c:102:     kd_node* pivot = arr[high];
	.loc 1 102 25 is_stmt 0 view .LVU329
	leaq	(%r9,%r13,8), %rbx	#, _45
# src/kdtree.c:104:     int i = (low - 1);
	.loc 1 104 9 view .LVU330
	movl	%r12d, %esi	# i, i
# src/kdtree.c:102:     kd_node* pivot = arr[high];
	.loc 1 102 14 view .LVU331
	movq	(%rbx), %rax	# *_45, pivot
.LVL86:
	.loc 1 104 5 is_stmt 1 view .LVU332
	.loc 1 105 5 view .LVU333
.LBB172:
	.loc 1 105 10 view .LVU334
	.loc 1 105 25 discriminator 1 view .LVU335
	cmpl	%r11d, %r10d	# right, left
	jge	.L75	#,
.L84:
	.loc 1 105 25 is_stmt 0 discriminator 1 view .LVU336
	movl	%r11d, %edx	# right, tmp166
.LBB167:
.LBB164:
# src/kdtree.c:79:     FLOAT_TYPE res = a->data[var] - b->data[var];
	.loc 1 79 44 view .LVU337
	movq	8(%rax), %rax	# pivot_46->data, pivot_46->data
.LVL87:
	.loc 1 79 44 view .LVU338
	movslq	%r10d, %rcx	# left, _89
	subl	%r10d, %edx	# left, tmp166
	addq	%rcx, %rdx	# _89, tmp167
	vmovsd	(%rax,%rdi), %xmm1	# *_59, _60
	leaq	(%r9,%rcx,8), %rax	#, ivtmp.139
	leaq	(%r9,%rdx,8), %r8	#, _150
.LVL88:
	.p2align 4,,10
	.p2align 3
.L69:
	.loc 1 79 44 view .LVU339
.LBE164:
.LBE167:
	.loc 1 106 9 is_stmt 1 view .LVU340
# src/kdtree.c:106:         if (!cmpKDnodes(arr[j],pivot,split_var)) {
	.loc 1 106 14 is_stmt 0 view .LVU341
	movq	(%rax), %rdx	# MEM[(struct kd_node * *)_50], _127
.LVL89:
.LBB168:
.LBI161:
	.loc 1 77 5 is_stmt 1 view .LVU342
.LBB165:
	.loc 1 79 5 view .LVU343
	.loc 1 80 5 view .LVU344
# src/kdtree.c:79:     FLOAT_TYPE res = a->data[var] - b->data[var];
	.loc 1 79 29 is_stmt 0 view .LVU345
	movq	8(%rdx), %rcx	# _127->data, _127->data
# src/kdtree.c:79:     FLOAT_TYPE res = a->data[var] - b->data[var];
	.loc 1 79 16 view .LVU346
	vmovsd	(%rcx,%rdi), %xmm0	# *_125, *_125
	vsubsd	%xmm1, %xmm0, %xmm0	# _60, *_125, res
.LBE165:
.LBE168:
# src/kdtree.c:106:         if (!cmpKDnodes(arr[j],pivot,split_var)) {
	.loc 1 106 12 discriminator 1 view .LVU347
	vcomisd	%xmm2, %xmm0	# tmp181, res
	ja	.L68	#,
# src/kdtree.c:107:             i++;
	.loc 1 107 14 view .LVU348
	addl	$1, %esi	#, i
.LVL90:
	.loc 1 107 13 is_stmt 1 view .LVU349
	.loc 1 108 13 view .LVU350
# src/kdtree.c:108:             swap_kd_node_ptrs(arr + i, arr + j);
	.loc 1 108 35 is_stmt 0 view .LVU351
	movslq	%esi, %rcx	# i, i
# src/kdtree.c:108:             swap_kd_node_ptrs(arr + i, arr + j);
	.loc 1 108 13 view .LVU352
	leaq	(%r9,%rcx,8), %rcx	#, _119
.LVL91:
.LBB169:
.LBI169:
	.loc 1 37 6 is_stmt 1 view .LVU353
.LBB170:
	.loc 1 38 5 view .LVU354
	.loc 1 39 5 view .LVU355
# src/kdtree.c:39:     tmp = *x;
	.loc 1 39 9 is_stmt 0 view .LVU356
	movq	(%rcx), %r14	# *_119, tmp
.LVL92:
	.loc 1 40 5 is_stmt 1 view .LVU357
# src/kdtree.c:40:     *x = *y;
	.loc 1 40 8 is_stmt 0 view .LVU358
	movq	%rdx, (%rcx)	# _127, *_119
	.loc 1 41 5 is_stmt 1 view .LVU359
# src/kdtree.c:41:     *y = tmp;
	.loc 1 41 8 is_stmt 0 view .LVU360
	movq	%r14, (%rax)	# tmp, MEM[(struct kd_node * *)_50]
.LVL93:
.L68:
	.loc 1 41 8 view .LVU361
.LBE170:
.LBE169:
	.loc 1 105 39 is_stmt 1 discriminator 2 view .LVU362
	.loc 1 105 25 discriminator 1 view .LVU363
	addq	$8, %rax	#, ivtmp.139
	cmpq	%r8, %rax	# _150, ivtmp.139
	jne	.L69	#,
.LVL94:
	.loc 1 105 25 is_stmt 0 discriminator 1 view .LVU364
.LBE172:
.LBB173:
.LBB174:
# src/kdtree.c:40:     *x = *y;
	.loc 1 40 10 view .LVU365
	movq	(%rbx), %rax	# *_45, pivot
.LBE174:
.LBE173:
# src/kdtree.c:112:     return (i + 1);
	.loc 1 112 15 view .LVU366
	leal	1(%rsi), %r8d	#, <retval>
.LVL95:
.L67:
	.loc 1 111 5 is_stmt 1 view .LVU367
# src/kdtree.c:111:     swap_kd_node_ptrs(arr + i + 1, arr + high);
	.loc 1 111 31 is_stmt 0 view .LVU368
	movslq	%esi, %rcx	# i, i
# src/kdtree.c:111:     swap_kd_node_ptrs(arr + i + 1, arr + high);
	.loc 1 111 5 view .LVU369
	leaq	8(%r9,%rcx,8), %rdx	#, _73
.LVL96:
.LBB176:
.LBI173:
	.loc 1 37 6 is_stmt 1 view .LVU370
.LBB175:
	.loc 1 38 5 view .LVU371
	.loc 1 39 5 view .LVU372
# src/kdtree.c:39:     tmp = *x;
	.loc 1 39 9 is_stmt 0 view .LVU373
	movq	(%rdx), %r14	# *_73, tmp
.LVL97:
	.loc 1 40 5 is_stmt 1 view .LVU374
# src/kdtree.c:40:     *x = *y;
	.loc 1 40 8 is_stmt 0 view .LVU375
	movq	%rax, (%rdx)	# pivot, *_73
	.loc 1 41 5 is_stmt 1 view .LVU376
# src/kdtree.c:41:     *y = tmp;
	.loc 1 41 8 is_stmt 0 view .LVU377
	movq	%r14, (%rbx)	# tmp, *_45
.LVL98:
	.loc 1 41 8 view .LVU378
.LBE175:
.LBE176:
	.loc 1 112 5 is_stmt 1 view .LVU379
	.loc 1 112 5 is_stmt 0 view .LVU380
.LBE179:
.LBE182:
	.loc 1 147 9 is_stmt 1 view .LVU381
# src/kdtree.c:147:         if (pivotIndex == k)
	.loc 1 147 12 is_stmt 0 view .LVU382
	cmpl	%r8d, %ebp	# <retval>, k
	je	.L61	#,
	.loc 1 153 14 is_stmt 1 view .LVU383
# src/kdtree.c:153:         else if (pivotIndex > k)
	.loc 1 153 17 is_stmt 0 view .LVU384
	jl	.L70	#,
	.loc 1 158 13 is_stmt 1 view .LVU385
# src/kdtree.c:158:             left = pivotIndex + 1;
	.loc 1 158 18 is_stmt 0 view .LVU386
	leal	2(%rsi), %r10d	#, left
.LVL99:
	.loc 1 158 18 view .LVU387
.LBE185:
	.loc 1 139 17 is_stmt 1 view .LVU388
	cmpl	%r11d, %r10d	# right, left
	jg	.L72	#,
	leal	1(%rsi), %r12d	#, i
.LVL100:
.LBB186:
	.loc 1 143 9 view .LVU389
.LBB183:
	.loc 1 100 5 view .LVU390
.LBB180:
	.loc 1 102 5 view .LVU391
# src/kdtree.c:102:     kd_node* pivot = arr[high];
	.loc 1 102 25 is_stmt 0 view .LVU392
	leaq	(%r9,%r13,8), %rbx	#, _45
# src/kdtree.c:102:     kd_node* pivot = arr[high];
	.loc 1 102 14 view .LVU393
	movq	(%rbx), %rax	# *_45, pivot
.LVL101:
	.loc 1 104 5 is_stmt 1 view .LVU394
# src/kdtree.c:104:     int i = (low - 1);
	.loc 1 104 9 is_stmt 0 view .LVU395
	movl	%r12d, %esi	# i, i
.LVL102:
	.loc 1 105 5 is_stmt 1 view .LVU396
.LBB177:
	.loc 1 105 10 view .LVU397
	.loc 1 105 25 discriminator 1 view .LVU398
	cmpl	%r11d, %r10d	# right, left
	jl	.L84	#,
.LVL103:
.L75:
	.loc 1 105 25 is_stmt 0 discriminator 1 view .LVU399
	movl	%r10d, %r8d	# left, <retval>
	jmp	.L67	#
.LVL104:
	.p2align 4,,10
	.p2align 3
.L70:
	.loc 1 105 25 discriminator 1 view .LVU400
.LBE177:
.LBE180:
.LBE183:
.LBE186:
	.loc 1 139 17 is_stmt 1 view .LVU401
	cmpl	%esi, %r10d	# i, left
	jg	.L72	#,
	.loc 1 139 17 is_stmt 0 view .LVU402
	movl	%esi, %r11d	# i, right
	movq	%rcx, %r13	# i, right
	jmp	.L73	#
.LVL105:
.L72:
# src/kdtree.c:161:     return -1;
	.loc 1 161 12 view .LVU403
	movl	$-1, %r8d	#, <retval>
.L61:
# src/kdtree.c:162: }
	.loc 1 162 1 view .LVU404
	popq	%rbx	#
	.cfi_def_cfa_offset 40
	movl	%r8d, %eax	# <retval>,
	popq	%rbp	#
	.cfi_def_cfa_offset 32
	popq	%r12	#
	.cfi_def_cfa_offset 24
	popq	%r13	#
	.cfi_def_cfa_offset 16
	popq	%r14	#
	.cfi_def_cfa_offset 8
	ret	
.LVL106:
.L63:
	.cfi_restore 3
	.cfi_restore 6
	.cfi_restore 12
	.cfi_restore 13
	.cfi_restore 14
.LBB187:
.LBI187:
	.loc 1 116 5 is_stmt 1 view .LVU405
.LBB188:
	.loc 1 131 9 view .LVU406
# src/kdtree.c:131:         if(cmpKDnodes(a[left],a[right],split_var)) {swap_kd_node_ptrs(a + left, a + right);}
	.loc 1 131 24 is_stmt 0 view .LVU407
	cltq
.LVL107:
# src/kdtree.c:131:         if(cmpKDnodes(a[left],a[right],split_var)) {swap_kd_node_ptrs(a + left, a + right);}
	.loc 1 131 32 view .LVU408
	movslq	%edx, %rdx	# right, right
.LVL108:
# src/kdtree.c:131:         if(cmpKDnodes(a[left],a[right],split_var)) {swap_kd_node_ptrs(a + left, a + right);}
	.loc 1 131 11 discriminator 1 view .LVU409
	vxorpd	%xmm1, %xmm1, %xmm1	# tmp156
.LBB189:
.LBB190:
# src/kdtree.c:79:     FLOAT_TYPE res = a->data[var] - b->data[var];
	.loc 1 79 29 view .LVU410
	movslq	%ecx, %rcx	# split_var, split_var
	.loc 1 79 29 view .LVU411
.LBE190:
.LBE189:
# src/kdtree.c:131:         if(cmpKDnodes(a[left],a[right],split_var)) {swap_kd_node_ptrs(a + left, a + right);}
	.loc 1 131 24 view .LVU412
	leaq	(%r9,%rax,8), %rsi	#, _31
.LVL109:
# src/kdtree.c:131:         if(cmpKDnodes(a[left],a[right],split_var)) {swap_kd_node_ptrs(a + left, a + right);}
	.loc 1 131 32 view .LVU413
	leaq	(%rdi,%rdx,8), %rdx	#, _27
# src/kdtree.c:131:         if(cmpKDnodes(a[left],a[right],split_var)) {swap_kd_node_ptrs(a + left, a + right);}
	.loc 1 131 12 view .LVU414
	movq	(%rdx), %rdi	# *_27, _28
.LVL110:
	.loc 1 131 12 view .LVU415
	movq	(%rsi), %rax	# *_31, _32
.LVL111:
.LBB192:
.LBI189:
	.loc 1 77 5 is_stmt 1 view .LVU416
.LBB191:
	.loc 1 79 5 view .LVU417
	.loc 1 80 5 view .LVU418
# src/kdtree.c:79:     FLOAT_TYPE res = a->data[var] - b->data[var];
	.loc 1 79 44 is_stmt 0 view .LVU419
	movq	8(%rdi), %r8	# _28->data, _28->data
.LVL112:
# src/kdtree.c:79:     FLOAT_TYPE res = a->data[var] - b->data[var];
	.loc 1 79 29 view .LVU420
	movq	8(%rax), %r9	# _32->data, _32->data
.LVL113:
# src/kdtree.c:79:     FLOAT_TYPE res = a->data[var] - b->data[var];
	.loc 1 79 16 view .LVU421
	vmovsd	(%r9,%rcx,8), %xmm0	# *_36, *_36
	vsubsd	(%r8,%rcx,8), %xmm0, %xmm0	# *_39, *_36, res
.LBE191:
.LBE192:
# src/kdtree.c:131:         if(cmpKDnodes(a[left],a[right],split_var)) {swap_kd_node_ptrs(a + left, a + right);}
	.loc 1 131 11 discriminator 1 view .LVU422
	vcomisd	%xmm1, %xmm0	# tmp156, res
	jbe	.L65	#,
	.loc 1 131 53 is_stmt 1 discriminator 1 view .LVU423
.LVL114:
.LBB193:
.LBI193:
	.loc 1 37 6 view .LVU424
.LBB194:
	.loc 1 38 5 view .LVU425
	.loc 1 39 5 view .LVU426
	.loc 1 40 5 view .LVU427
# src/kdtree.c:40:     *x = *y;
	.loc 1 40 8 is_stmt 0 view .LVU428
	movq	%rdi, (%rsi)	# _28, *_31
	.loc 1 41 5 is_stmt 1 view .LVU429
# src/kdtree.c:41:     *y = tmp;
	.loc 1 41 8 is_stmt 0 view .LVU430
	movq	%rax, (%rdx)	# _32, *_27
.LVL115:
.L65:
	.loc 1 41 8 view .LVU431
.LBE194:
.LBE193:
.LBE188:
.LBE187:
# src/kdtree.c:129:     if(left == right) return left;
	.loc 1 129 30 discriminator 1 view .LVU432
	movl	%r11d, %r8d	# right, <retval>
.LVL116:
.L80:
# src/kdtree.c:162: }
	.loc 1 162 1 view .LVU433
	movl	%r8d, %eax	# <retval>,
	ret	
.LVL117:
.L83:
# src/kdtree.c:161:     return -1;
	.loc 1 161 12 view .LVU434
	movl	$-1, %r8d	#, <retval>
	jmp	.L80	#
	.cfi_endproc
.LFE62:
	.size	medianOfNodes, .-medianOfNodes
	.p2align 4
	.globl	make_tree
	.type	make_tree, @function
make_tree:
.LVL118:
.LFB63:
	.loc 1 165 1 is_stmt 1 view -0
	.cfi_startproc
	.loc 1 165 1 is_stmt 0 view .LVU436
	endbr64	
	.loc 1 166 5 is_stmt 1 view .LVU437
.LVL119:
	.loc 1 167 5 view .LVU438
# src/kdtree.c:165: {
	.loc 1 165 1 is_stmt 0 view .LVU439
	pushq	%r15	#
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
# src/kdtree.c:167:     int split_var = level % data_dims; 
	.loc 1 167 27 view .LVU440
	movl	%r8d, %eax	# level, tmp106
# src/kdtree.c:165: {
	.loc 1 165 1 view .LVU441
	movslq	%edx, %r15	# tmp114,
# src/kdtree.c:167:     int split_var = level % data_dims; 
	.loc 1 167 27 view .LVU442
	xorl	%edx, %edx	# tmp105
.LVL120:
# src/kdtree.c:165: {
	.loc 1 165 1 view .LVU443
	pushq	%r14	#
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	pushq	%r13	#
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	pushq	%r12	#
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	pushq	%rbp	#
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	pushq	%rbx	#
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	subq	$24, %rsp	#,
	.cfi_def_cfa_offset 80
# src/kdtree.c:167:     int split_var = level % data_dims; 
	.loc 1 167 27 view .LVU444
	divl	data_dims(%rip)	# data_dims
.LVL121:
	.loc 1 170 5 is_stmt 1 view .LVU445
	.loc 1 171 5 view .LVU446
# src/kdtree.c:171:     if ((end - start) < 0) return 0;
	.loc 1 171 8 is_stmt 0 view .LVU447
	cmpl	%esi, %r15d	# start, end
	jl	.L89	#,
	movq	%rcx, %r13	# tmp115, parent
	movl	%r8d, %ebp	# tmp116, level
	movl	%edx, %r14d	# tmp105, tmp105
	.loc 1 172 5 is_stmt 1 view .LVU448
# src/kdtree.c:172:     if (end  == start) {
	.loc 1 172 8 is_stmt 0 view .LVU449
	je	.L91	#,
	.loc 1 181 5 is_stmt 1 view .LVU450
# src/kdtree.c:181:     median_idx = medianOfNodes(t, start, end, split_var);
	.loc 1 181 18 is_stmt 0 view .LVU451
	movl	%edx, %ecx	# tmp105,
.LVL122:
	.loc 1 181 18 view .LVU452
	movl	%r15d, %edx	# end,
.LVL123:
	.loc 1 181 18 view .LVU453
	movl	%esi, 8(%rsp)	# start, %sfp
	movq	%rdi, (%rsp)	# t, %sfp
	call	medianOfNodes	#
.LVL124:
	.loc 1 181 18 view .LVU454
	movl	%eax, %r12d	# tmp117, median_idx
.LVL125:
	.loc 1 183 5 is_stmt 1 view .LVU455
# src/kdtree.c:183:     if(median_idx > -1){
	.loc 1 183 7 is_stmt 0 view .LVU456
	testl	%eax, %eax	# median_idx
	js	.L89	#,
	.loc 1 184 9 is_stmt 1 view .LVU457
# src/kdtree.c:184:         n = t[median_idx];
	.loc 1 184 11 is_stmt 0 view .LVU458
	movq	(%rsp), %rdi	# %sfp, t
# src/kdtree.c:184:         n = t[median_idx];
	.loc 1 184 14 view .LVU459
	cltq
.LVL126:
# src/kdtree.c:190: 				n->lch  = make_tree(t, start, median_idx - 1, n, level + 1);
	.loc 1 190 15 view .LVU460
	movl	8(%rsp), %esi	# %sfp, start
	leal	1(%rbp), %r8d	#, _11
	leal	-1(%r12), %edx	#, tmp109
	movl	%r8d, 12(%rsp)	# _11, %sfp
# src/kdtree.c:184:         n = t[median_idx];
	.loc 1 184 11 view .LVU461
	movq	(%rdi,%rax,8), %rbx	# *_10, <retval>
.LVL127:
	.loc 1 190 5 is_stmt 1 view .LVU462
# src/kdtree.c:190: 				n->lch  = make_tree(t, start, median_idx - 1, n, level + 1);
	.loc 1 190 15 is_stmt 0 view .LVU463
	movq	%rbx, %rcx	# <retval>,
	call	make_tree	#
.LVL128:
# src/kdtree.c:192: 				n->rch = make_tree(t, median_idx + 1, end, n, level + 1);
	.loc 1 192 14 view .LVU464
	movl	12(%rsp), %r8d	# %sfp, _11
	movq	%rbx, %rcx	# <retval>,
	movl	%r15d, %edx	# end,
# src/kdtree.c:190: 				n->lch  = make_tree(t, start, median_idx - 1, n, level + 1);
	.loc 1 190 13 discriminator 1 view .LVU465
	movq	%rax, 32(%rbx)	# tmp118, n_27->lch
	.loc 1 192 5 is_stmt 1 view .LVU466
# src/kdtree.c:192: 				n->rch = make_tree(t, median_idx + 1, end, n, level + 1);
	.loc 1 192 14 is_stmt 0 view .LVU467
	movq	(%rsp), %rdi	# %sfp, t
	leal	1(%r12), %esi	#, tmp110
	call	make_tree	#
.LVL129:
# src/kdtree.c:195:         n -> split_var = split_var;
	.loc 1 195 24 view .LVU468
	movl	%r14d, 4(%rbx)	# tmp105, n_27->split_var
# src/kdtree.c:192: 				n->rch = make_tree(t, median_idx + 1, end, n, level + 1);
	.loc 1 192 12 discriminator 1 view .LVU469
	movq	%rax, 40(%rbx)	# tmp119, n_27->rch
	.loc 1 195 9 is_stmt 1 view .LVU470
	.loc 1 196 9 view .LVU471
# src/kdtree.c:196:         n->parent = parent;
	.loc 1 196 19 is_stmt 0 view .LVU472
	movq	%r13, 24(%rbx)	# parent, n_27->parent
	.loc 1 197 9 is_stmt 1 view .LVU473
# src/kdtree.c:197:         n->level = level;
	.loc 1 197 18 is_stmt 0 view .LVU474
	movl	%ebp, (%rbx)	# level, n_27->level
.LVL130:
.L85:
# src/kdtree.c:210: }
	.loc 1 210 1 view .LVU475
	addq	$24, %rsp	#,
	.cfi_remember_state
	.cfi_def_cfa_offset 56
	movq	%rbx, %rax	# <retval>,
	popq	%rbx	#
	.cfi_def_cfa_offset 48
	popq	%rbp	#
	.cfi_def_cfa_offset 40
	popq	%r12	#
	.cfi_def_cfa_offset 32
	popq	%r13	#
	.cfi_def_cfa_offset 24
	popq	%r14	#
	.cfi_def_cfa_offset 16
	popq	%r15	#
	.cfi_def_cfa_offset 8
.LVL131:
	.loc 1 210 1 view .LVU476
	ret	
.LVL132:
	.p2align 4,,10
	.p2align 3
.L89:
	.cfi_restore_state
# src/kdtree.c:171:     if ((end - start) < 0) return 0;
	.loc 1 171 35 discriminator 1 view .LVU477
	xorl	%ebx, %ebx	# <retval>
	jmp	.L85	#
.LVL133:
	.p2align 4,,10
	.p2align 3
.L91:
	.loc 1 173 9 is_stmt 1 view .LVU478
# src/kdtree.c:173:         n = t[start];
	.loc 1 173 11 is_stmt 0 view .LVU479
	movq	(%rdi,%r15,8), %rbx	# *_7, <retval>
.LVL134:
	.loc 1 174 9 is_stmt 1 view .LVU480
# src/kdtree.c:174:         n -> split_var = split_var;
	.loc 1 174 24 is_stmt 0 view .LVU481
	movl	%edx, 4(%rbx)	# tmp105, n_36->split_var
	.loc 1 175 9 is_stmt 1 view .LVU482
# src/kdtree.c:175:         n->parent = parent;
	.loc 1 175 19 is_stmt 0 view .LVU483
	movq	%rcx, 24(%rbx)	# parent, n_36->parent
	.loc 1 176 9 is_stmt 1 view .LVU484
# src/kdtree.c:176:         n->level = level;
	.loc 1 176 18 is_stmt 0 view .LVU485
	movl	%r8d, (%rbx)	# level, n_36->level
	.loc 1 177 9 is_stmt 1 view .LVU486
# src/kdtree.c:177:         n -> lch = NULL;
	.loc 1 177 18 is_stmt 0 view .LVU487
	movq	$0, 32(%rbx)	#, n_36->lch
	.loc 1 178 9 is_stmt 1 view .LVU488
# src/kdtree.c:178:         n -> rch = NULL;
	.loc 1 178 18 is_stmt 0 view .LVU489
	movq	$0, 40(%rbx)	#, n_36->rch
	.loc 1 179 9 is_stmt 1 view .LVU490
# src/kdtree.c:179:         return n;
	.loc 1 179 16 is_stmt 0 view .LVU491
	jmp	.L85	#
	.cfi_endproc
.LFE63:
	.size	make_tree, .-make_tree
	.p2align 4
	.type	make_tree.constprop.2, @function
make_tree.constprop.2:
.LVL135:
.LFB71:
	.loc 1 164 10 is_stmt 1 view -0
	.cfi_startproc
	.loc 1 166 5 view .LVU493
	.loc 1 167 5 view .LVU494
# src/kdtree.c:164: kd_node* make_tree(kd_node** t, int start, int end, kd_node* parent, int level)
	.loc 1 164 10 is_stmt 0 view .LVU495
	pushq	%r15	#
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
# src/kdtree.c:167:     int split_var = level % data_dims; 
	.loc 1 167 27 view .LVU496
	movl	$2, %eax	#, tmp102
# src/kdtree.c:164: kd_node* make_tree(kd_node** t, int start, int end, kd_node* parent, int level)
	.loc 1 164 10 view .LVU497
	pushq	%r14	#
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	pushq	%r13	#
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	pushq	%r12	#
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	pushq	%rbp	#
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	movslq	%edx, %rbp	# tmp112,
# src/kdtree.c:167:     int split_var = level % data_dims; 
	.loc 1 167 27 view .LVU498
	xorl	%edx, %edx	# tmp103
.LVL136:
# src/kdtree.c:164: kd_node* make_tree(kd_node** t, int start, int end, kd_node* parent, int level)
	.loc 1 164 10 view .LVU499
	pushq	%rbx	#
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	subq	$24, %rsp	#,
	.cfi_def_cfa_offset 80
# src/kdtree.c:167:     int split_var = level % data_dims; 
	.loc 1 167 27 view .LVU500
	divl	data_dims(%rip)	# data_dims
.LVL137:
	.loc 1 170 5 is_stmt 1 view .LVU501
	.loc 1 171 5 view .LVU502
# src/kdtree.c:171:     if ((end - start) < 0) return 0;
	.loc 1 171 8 is_stmt 0 view .LVU503
	cmpl	%esi, %ebp	# start, end
	jl	.L96	#,
	movq	%rdi, %r12	# tmp110, t
	movq	%rcx, %r14	# tmp113, parent
	movl	%edx, %r15d	# tmp103, tmp103
	.loc 1 172 5 is_stmt 1 view .LVU504
# src/kdtree.c:172:     if (end  == start) {
	.loc 1 172 8 is_stmt 0 view .LVU505
	je	.L98	#,
# src/kdtree.c:181:     median_idx = medianOfNodes(t, start, end, split_var);
	.loc 1 181 18 view .LVU506
	movl	%edx, %ecx	# tmp103,
.LVL138:
	.loc 1 181 18 view .LVU507
	movl	%ebp, %edx	# end,
.LVL139:
	.loc 1 181 18 view .LVU508
	movl	%esi, %r13d	# tmp111, start
	.loc 1 181 5 is_stmt 1 view .LVU509
# src/kdtree.c:181:     median_idx = medianOfNodes(t, start, end, split_var);
	.loc 1 181 18 is_stmt 0 view .LVU510
	call	medianOfNodes	#
.LVL140:
	.loc 1 181 18 view .LVU511
	movl	%eax, %r9d	# tmp114, median_idx
.LVL141:
	.loc 1 183 5 is_stmt 1 view .LVU512
# src/kdtree.c:183:     if(median_idx > -1){
	.loc 1 183 7 is_stmt 0 view .LVU513
	testl	%eax, %eax	# median_idx
	js	.L96	#,
	.loc 1 184 9 is_stmt 1 view .LVU514
# src/kdtree.c:184:         n = t[median_idx];
	.loc 1 184 14 is_stmt 0 view .LVU515
	cltq
.LVL142:
# src/kdtree.c:190: 				n->lch  = make_tree(t, start, median_idx - 1, n, level + 1);
	.loc 1 190 15 view .LVU516
	leal	-1(%r9), %edx	#, tmp107
	movl	$3, %r8d	#,
	movl	%r13d, %esi	# start,
# src/kdtree.c:184:         n = t[median_idx];
	.loc 1 184 11 view .LVU517
	movq	(%r12,%rax,8), %rbx	# *_17, <retval>
.LVL143:
	.loc 1 190 5 is_stmt 1 view .LVU518
# src/kdtree.c:190: 				n->lch  = make_tree(t, start, median_idx - 1, n, level + 1);
	.loc 1 190 15 is_stmt 0 view .LVU519
	movq	%r12, %rdi	# t,
	movl	%r9d, 12(%rsp)	# median_idx, %sfp
	movq	%rbx, %rcx	# <retval>,
	call	make_tree	#
.LVL144:
# src/kdtree.c:192: 				n->rch = make_tree(t, median_idx + 1, end, n, level + 1);
	.loc 1 192 14 view .LVU520
	movl	12(%rsp), %r9d	# %sfp, median_idx
	movq	%rbx, %rcx	# <retval>,
	movl	%ebp, %edx	# end,
# src/kdtree.c:190: 				n->lch  = make_tree(t, start, median_idx - 1, n, level + 1);
	.loc 1 190 13 discriminator 1 view .LVU521
	movq	%rax, 32(%rbx)	# tmp115, n_18->lch
	.loc 1 192 5 is_stmt 1 view .LVU522
# src/kdtree.c:192: 				n->rch = make_tree(t, median_idx + 1, end, n, level + 1);
	.loc 1 192 14 is_stmt 0 view .LVU523
	movl	$3, %r8d	#,
	movq	%r12, %rdi	# t,
	leal	1(%r9), %esi	#, tmp108
	call	make_tree	#
.LVL145:
# src/kdtree.c:195:         n -> split_var = split_var;
	.loc 1 195 24 view .LVU524
	movl	%r15d, 4(%rbx)	# tmp103, n_18->split_var
# src/kdtree.c:192: 				n->rch = make_tree(t, median_idx + 1, end, n, level + 1);
	.loc 1 192 12 discriminator 1 view .LVU525
	movq	%rax, 40(%rbx)	# tmp116, n_18->rch
	.loc 1 195 9 is_stmt 1 view .LVU526
	.loc 1 196 9 view .LVU527
# src/kdtree.c:196:         n->parent = parent;
	.loc 1 196 19 is_stmt 0 view .LVU528
	movq	%r14, 24(%rbx)	# parent, n_18->parent
	.loc 1 197 9 is_stmt 1 view .LVU529
# src/kdtree.c:197:         n->level = level;
	.loc 1 197 18 is_stmt 0 view .LVU530
	movl	$2, (%rbx)	#, n_18->level
.LVL146:
.L92:
# src/kdtree.c:210: }
	.loc 1 210 1 view .LVU531
	addq	$24, %rsp	#,
	.cfi_remember_state
	.cfi_def_cfa_offset 56
	movq	%rbx, %rax	# <retval>,
	popq	%rbx	#
	.cfi_def_cfa_offset 48
	popq	%rbp	#
	.cfi_def_cfa_offset 40
.LVL147:
	.loc 1 210 1 view .LVU532
	popq	%r12	#
	.cfi_def_cfa_offset 32
	popq	%r13	#
	.cfi_def_cfa_offset 24
	popq	%r14	#
	.cfi_def_cfa_offset 16
	popq	%r15	#
	.cfi_def_cfa_offset 8
	ret	
.LVL148:
	.p2align 4,,10
	.p2align 3
.L96:
	.cfi_restore_state
# src/kdtree.c:171:     if ((end - start) < 0) return 0;
	.loc 1 171 35 discriminator 1 view .LVU533
	xorl	%ebx, %ebx	# <retval>
	jmp	.L92	#
.LVL149:
	.p2align 4,,10
	.p2align 3
.L98:
	.loc 1 173 9 is_stmt 1 view .LVU534
# src/kdtree.c:173:         n = t[start];
	.loc 1 173 11 is_stmt 0 view .LVU535
	movq	(%rdi,%rbp,8), %rbx	# *_11, <retval>
.LVL150:
	.loc 1 174 9 is_stmt 1 view .LVU536
# src/kdtree.c:174:         n -> split_var = split_var;
	.loc 1 174 24 is_stmt 0 view .LVU537
	movl	%edx, 4(%rbx)	# tmp103, n_12->split_var
	.loc 1 175 9 is_stmt 1 view .LVU538
# src/kdtree.c:175:         n->parent = parent;
	.loc 1 175 19 is_stmt 0 view .LVU539
	movq	%rcx, 24(%rbx)	# parent, n_12->parent
	.loc 1 176 9 is_stmt 1 view .LVU540
# src/kdtree.c:176:         n->level = level;
	.loc 1 176 18 is_stmt 0 view .LVU541
	movl	$2, (%rbx)	#, n_12->level
	.loc 1 177 9 is_stmt 1 view .LVU542
# src/kdtree.c:177:         n -> lch = NULL;
	.loc 1 177 18 is_stmt 0 view .LVU543
	movq	$0, 32(%rbx)	#, n_12->lch
	.loc 1 178 9 is_stmt 1 view .LVU544
# src/kdtree.c:178:         n -> rch = NULL;
	.loc 1 178 18 is_stmt 0 view .LVU545
	movq	$0, 40(%rbx)	#, n_12->rch
	.loc 1 179 9 is_stmt 1 view .LVU546
# src/kdtree.c:179:         return n;
	.loc 1 179 16 is_stmt 0 view .LVU547
	jmp	.L92	#
	.cfi_endproc
.LFE71:
	.size	make_tree.constprop.2, .-make_tree.constprop.2
	.p2align 4
	.type	make_tree.constprop.1, @function
make_tree.constprop.1:
.LVL151:
.LFB72:
	.loc 1 164 10 is_stmt 1 view -0
	.cfi_startproc
	.loc 1 166 5 view .LVU549
	.loc 1 167 5 view .LVU550
# src/kdtree.c:164: kd_node* make_tree(kd_node** t, int start, int end, kd_node* parent, int level)
	.loc 1 164 10 is_stmt 0 view .LVU551
	pushq	%r15	#
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
# src/kdtree.c:167:     int split_var = level % data_dims; 
	.loc 1 167 27 view .LVU552
	movl	$1, %eax	#, tmp102
# src/kdtree.c:164: kd_node* make_tree(kd_node** t, int start, int end, kd_node* parent, int level)
	.loc 1 164 10 view .LVU553
	pushq	%r14	#
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	pushq	%r13	#
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	pushq	%r12	#
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	pushq	%rbp	#
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	movslq	%edx, %rbp	# tmp112,
# src/kdtree.c:167:     int split_var = level % data_dims; 
	.loc 1 167 27 view .LVU554
	xorl	%edx, %edx	# tmp103
.LVL152:
# src/kdtree.c:164: kd_node* make_tree(kd_node** t, int start, int end, kd_node* parent, int level)
	.loc 1 164 10 view .LVU555
	pushq	%rbx	#
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	subq	$24, %rsp	#,
	.cfi_def_cfa_offset 80
# src/kdtree.c:167:     int split_var = level % data_dims; 
	.loc 1 167 27 view .LVU556
	divl	data_dims(%rip)	# data_dims
.LVL153:
	.loc 1 170 5 is_stmt 1 view .LVU557
	.loc 1 171 5 view .LVU558
# src/kdtree.c:171:     if ((end - start) < 0) return 0;
	.loc 1 171 8 is_stmt 0 view .LVU559
	cmpl	%esi, %ebp	# start, end
	jl	.L103	#,
	movq	%rdi, %r12	# tmp110, t
	movq	%rcx, %r14	# tmp113, parent
	movl	%edx, %r15d	# tmp103, tmp103
	.loc 1 172 5 is_stmt 1 view .LVU560
# src/kdtree.c:172:     if (end  == start) {
	.loc 1 172 8 is_stmt 0 view .LVU561
	je	.L105	#,
# src/kdtree.c:181:     median_idx = medianOfNodes(t, start, end, split_var);
	.loc 1 181 18 view .LVU562
	movl	%edx, %ecx	# tmp103,
.LVL154:
	.loc 1 181 18 view .LVU563
	movl	%ebp, %edx	# end,
.LVL155:
	.loc 1 181 18 view .LVU564
	movl	%esi, %r13d	# tmp111, start
	.loc 1 181 5 is_stmt 1 view .LVU565
# src/kdtree.c:181:     median_idx = medianOfNodes(t, start, end, split_var);
	.loc 1 181 18 is_stmt 0 view .LVU566
	call	medianOfNodes	#
.LVL156:
	.loc 1 181 18 view .LVU567
	movl	%eax, %r8d	# tmp114, median_idx
.LVL157:
	.loc 1 183 5 is_stmt 1 view .LVU568
# src/kdtree.c:183:     if(median_idx > -1){
	.loc 1 183 7 is_stmt 0 view .LVU569
	testl	%eax, %eax	# median_idx
	js	.L103	#,
	.loc 1 184 9 is_stmt 1 view .LVU570
# src/kdtree.c:184:         n = t[median_idx];
	.loc 1 184 14 is_stmt 0 view .LVU571
	cltq
.LVL158:
# src/kdtree.c:190: 				n->lch  = make_tree(t, start, median_idx - 1, n, level + 1);
	.loc 1 190 15 view .LVU572
	leal	-1(%r8), %edx	#, tmp107
	movl	%r13d, %esi	# start,
	movq	%r12, %rdi	# t,
# src/kdtree.c:184:         n = t[median_idx];
	.loc 1 184 11 view .LVU573
	movq	(%r12,%rax,8), %rbx	# *_17, <retval>
.LVL159:
	.loc 1 190 5 is_stmt 1 view .LVU574
# src/kdtree.c:190: 				n->lch  = make_tree(t, start, median_idx - 1, n, level + 1);
	.loc 1 190 15 is_stmt 0 view .LVU575
	movl	%r8d, 12(%rsp)	# median_idx, %sfp
	movq	%rbx, %rcx	# <retval>,
	call	make_tree.constprop.2	#
.LVL160:
# src/kdtree.c:192: 				n->rch = make_tree(t, median_idx + 1, end, n, level + 1);
	.loc 1 192 14 view .LVU576
	movl	12(%rsp), %r8d	# %sfp, median_idx
	movq	%rbx, %rcx	# <retval>,
	movl	%ebp, %edx	# end,
# src/kdtree.c:190: 				n->lch  = make_tree(t, start, median_idx - 1, n, level + 1);
	.loc 1 190 13 discriminator 1 view .LVU577
	movq	%rax, 32(%rbx)	# tmp115, n_18->lch
	.loc 1 192 5 is_stmt 1 view .LVU578
# src/kdtree.c:192: 				n->rch = make_tree(t, median_idx + 1, end, n, level + 1);
	.loc 1 192 14 is_stmt 0 view .LVU579
	movq	%r12, %rdi	# t,
	leal	1(%r8), %esi	#, tmp108
	call	make_tree.constprop.2	#
.LVL161:
# src/kdtree.c:195:         n -> split_var = split_var;
	.loc 1 195 24 view .LVU580
	movl	%r15d, 4(%rbx)	# tmp103, n_18->split_var
# src/kdtree.c:192: 				n->rch = make_tree(t, median_idx + 1, end, n, level + 1);
	.loc 1 192 12 discriminator 1 view .LVU581
	movq	%rax, 40(%rbx)	# tmp116, n_18->rch
	.loc 1 195 9 is_stmt 1 view .LVU582
	.loc 1 196 9 view .LVU583
# src/kdtree.c:196:         n->parent = parent;
	.loc 1 196 19 is_stmt 0 view .LVU584
	movq	%r14, 24(%rbx)	# parent, n_18->parent
	.loc 1 197 9 is_stmt 1 view .LVU585
# src/kdtree.c:197:         n->level = level;
	.loc 1 197 18 is_stmt 0 view .LVU586
	movl	$1, (%rbx)	#, n_18->level
.LVL162:
.L99:
# src/kdtree.c:210: }
	.loc 1 210 1 view .LVU587
	addq	$24, %rsp	#,
	.cfi_remember_state
	.cfi_def_cfa_offset 56
	movq	%rbx, %rax	# <retval>,
	popq	%rbx	#
	.cfi_def_cfa_offset 48
	popq	%rbp	#
	.cfi_def_cfa_offset 40
.LVL163:
	.loc 1 210 1 view .LVU588
	popq	%r12	#
	.cfi_def_cfa_offset 32
	popq	%r13	#
	.cfi_def_cfa_offset 24
	popq	%r14	#
	.cfi_def_cfa_offset 16
	popq	%r15	#
	.cfi_def_cfa_offset 8
	ret	
.LVL164:
	.p2align 4,,10
	.p2align 3
.L103:
	.cfi_restore_state
# src/kdtree.c:171:     if ((end - start) < 0) return 0;
	.loc 1 171 35 discriminator 1 view .LVU589
	xorl	%ebx, %ebx	# <retval>
	jmp	.L99	#
.LVL165:
	.p2align 4,,10
	.p2align 3
.L105:
	.loc 1 173 9 is_stmt 1 view .LVU590
# src/kdtree.c:173:         n = t[start];
	.loc 1 173 11 is_stmt 0 view .LVU591
	movq	(%rdi,%rbp,8), %rbx	# *_11, <retval>
.LVL166:
	.loc 1 174 9 is_stmt 1 view .LVU592
# src/kdtree.c:174:         n -> split_var = split_var;
	.loc 1 174 24 is_stmt 0 view .LVU593
	movl	%edx, 4(%rbx)	# tmp103, n_12->split_var
	.loc 1 175 9 is_stmt 1 view .LVU594
# src/kdtree.c:175:         n->parent = parent;
	.loc 1 175 19 is_stmt 0 view .LVU595
	movq	%rcx, 24(%rbx)	# parent, n_12->parent
	.loc 1 176 9 is_stmt 1 view .LVU596
# src/kdtree.c:176:         n->level = level;
	.loc 1 176 18 is_stmt 0 view .LVU597
	movl	$1, (%rbx)	#, n_12->level
	.loc 1 177 9 is_stmt 1 view .LVU598
# src/kdtree.c:177:         n -> lch = NULL;
	.loc 1 177 18 is_stmt 0 view .LVU599
	movq	$0, 32(%rbx)	#, n_12->lch
	.loc 1 178 9 is_stmt 1 view .LVU600
# src/kdtree.c:178:         n -> rch = NULL;
	.loc 1 178 18 is_stmt 0 view .LVU601
	movq	$0, 40(%rbx)	#, n_12->rch
	.loc 1 179 9 is_stmt 1 view .LVU602
# src/kdtree.c:179:         return n;
	.loc 1 179 16 is_stmt 0 view .LVU603
	jmp	.L99	#
	.cfi_endproc
.LFE72:
	.size	make_tree.constprop.1, .-make_tree.constprop.1
	.p2align 4
	.globl	KNN_sub_tree_search
	.type	KNN_sub_tree_search, @function
KNN_sub_tree_search:
.LVL167:
.LFB66:
	.loc 1 223 1 is_stmt 1 view -0
	.cfi_startproc
	.loc 1 223 1 is_stmt 0 view .LVU605
	endbr64	
	leaq	8(%rsp), %r10	#,
	.cfi_def_cfa 10, 0
	andq	$-32, %rsp	#,
	pushq	-8(%r10)	#
	pushq	%rbp	#
	movq	%rsp, %rbp	#,
	.cfi_escape 0x10,0x6,0x2,0x76,0
	pushq	%r13	#
	.cfi_escape 0x10,0xd,0x2,0x76,0x78
	movq	%rdx, %r13	# tmp200, H
	pushq	%r12	#
	.cfi_escape 0x10,0xc,0x2,0x76,0x70
	movq	%rsi, %r12	# tmp199, kdtree_root
	pushq	%r10	#
	.cfi_escape 0xf,0x3,0x76,0x68,0x6
	pushq	%rbx	#
	.cfi_escape 0x10,0x3,0x2,0x76,0x60
	movq	%rdi, %rbx	# tmp198, point
	subq	$48, %rsp	#,
.LVL168:
.L107:
	.loc 1 225 5 is_stmt 1 view .LVU606
.LBB195:
.LBB196:
.LBB197:
# src/kdtree.c:21:     for(unsigned int i = 0; i<data_dims; ++i){
	.loc 1 21 30 is_stmt 0 discriminator 1 view .LVU607
	movl	data_dims(%rip), %ecx	# data_dims, data_dims.0_51
.LBE197:
.LBE196:
.LBE195:
# src/kdtree.c:225:     FLOAT_TYPE current_distance = euclidean_distance(point, kdtree_root -> data);
	.loc 1 225 35 view .LVU608
	movq	8(%r12), %rdx	# kdtree_root_18->data, _1
.LVL169:
.LBB210:
.LBI195:
	.loc 1 19 12 is_stmt 1 view .LVU609
.LBB206:
	.loc 1 20 5 view .LVU610
	.loc 1 21 5 view .LVU611
.LBB201:
	.loc 1 21 9 view .LVU612
	.loc 1 21 30 discriminator 1 view .LVU613
	testl	%ecx, %ecx	# data_dims.0_51
	je	.L124	#,
.L143:
	.loc 1 21 30 is_stmt 0 discriminator 1 view .LVU614
	leal	-1(%rcx), %eax	#, tmp157
	cmpl	$2, %eax	#, tmp157
	jbe	.L125	#,
	movl	%ecx, %esi	# data_dims.0_51, bnd.156
	xorl	%eax, %eax	# ivtmp.188
.LBE201:
# src/kdtree.c:20:     FLOAT_TYPE d = 0;
	.loc 1 20 16 view .LVU615
	vxorpd	%xmm0, %xmm0, %xmm0	# d
	shrl	$2, %esi	#,
	salq	$5, %rsi	#, _40
.LVL170:
	.p2align 4,,10
	.p2align 3
.L110:
.LBB202:
.LBB198:
	.loc 1 22 9 is_stmt 1 view .LVU616
# src/kdtree.c:22:         FLOAT_TYPE dd = (p1[i] - p2[i]);
	.loc 1 22 20 is_stmt 0 view .LVU617
	vmovupd	(%rbx,%rax), %ymm4	# MEM <vector(4) double> [(double *)point_19(D) + ivtmp.188_43 * 1], tmp203
	vsubpd	(%rdx,%rax), %ymm4, %ymm1	# MEM <vector(4) double> [(double *)_1 + ivtmp.188_43 * 1], tmp203, vect_dd_42.165
.LVL171:
	.loc 1 23 9 is_stmt 1 view .LVU618
	addq	$32, %rax	#, ivtmp.188
# src/kdtree.c:23:         d += dd*dd;
	.loc 1 23 16 is_stmt 0 view .LVU619
	vmulpd	%ymm1, %ymm1, %ymm1	# vect_dd_42.165, vect_dd_42.165, vect__43.166
	vaddsd	%xmm1, %xmm0, %xmm0	# stmp_d_45.167, d, stmp_d_45.167
.LVL172:
	.loc 1 23 16 view .LVU620
	vunpckhpd	%xmm1, %xmm1, %xmm2	# tmp162, stmp_d_45.167
	vextractf128	$0x1, %ymm1, %xmm1	# vect__43.166, tmp164
	vaddsd	%xmm0, %xmm2, %xmm2	# stmp_d_45.167, stmp_d_45.167, stmp_d_45.167
# src/kdtree.c:23:         d += dd*dd;
	.loc 1 23 11 view .LVU621
	vaddsd	%xmm2, %xmm1, %xmm0	# stmp_d_45.167, stmp_d_45.167, stmp_d_45.167
	vunpckhpd	%xmm1, %xmm1, %xmm1	# tmp164, stmp_d_45.167
	vaddsd	%xmm1, %xmm0, %xmm0	# stmp_d_45.167, stmp_d_45.167, d
.LVL173:
	.loc 1 23 11 view .LVU622
.LBE198:
	.loc 1 21 42 is_stmt 1 discriminator 3 view .LVU623
	.loc 1 21 30 discriminator 1 view .LVU624
	cmpq	%rax, %rsi	# ivtmp.188, _40
	jne	.L110	#,
	.loc 1 21 30 is_stmt 0 discriminator 1 view .LVU625
	testb	$3, %cl	#, data_dims.0_51
	je	.L136	#,
	movl	%ecx, %eax	# data_dims.0_51, niters_vector_mult_vf.157
	andl	$-4, %eax	#,
	vzeroupper
.LVL174:
.L109:
	.loc 1 21 30 discriminator 1 view .LVU626
	subl	%eax, %ecx	# niters_vector_mult_vf.157, niters.168
	cmpl	$1, %ecx	#, niters.168
	je	.L112	#,
	movl	%eax, %esi	# niters_vector_mult_vf.157, niters_vector_mult_vf.157
.LVL175:
.LBB199:
	.loc 1 22 9 is_stmt 1 view .LVU627
# src/kdtree.c:22:         FLOAT_TYPE dd = (p1[i] - p2[i]);
	.loc 1 22 20 is_stmt 0 view .LVU628
	vmovupd	(%rbx,%rsi,8), %xmm5	# MEM <vector(2) double> [(double *)vectp_point.173_83], tmp205
	vsubpd	(%rdx,%rsi,8), %xmm5, %xmm1	# MEM <vector(2) double> [(double *)vectp.176_77], tmp205, vect_dd_136.178
.LVL176:
	.loc 1 23 9 is_stmt 1 view .LVU629
# src/kdtree.c:23:         d += dd*dd;
	.loc 1 23 16 is_stmt 0 view .LVU630
	vmulpd	%xmm1, %xmm1, %xmm1	# vect_dd_136.178, vect_dd_136.178, vect__135.179
	vaddsd	%xmm0, %xmm1, %xmm0	# d, stmp_d_134.180, stmp_d_134.180
.LVL177:
# src/kdtree.c:23:         d += dd*dd;
	.loc 1 23 11 view .LVU631
	vunpckhpd	%xmm1, %xmm1, %xmm1	# vect__135.179, stmp_d_134.180
	vaddsd	%xmm0, %xmm1, %xmm0	# stmp_d_134.180, stmp_d_134.180, d
.LVL178:
	.loc 1 23 11 view .LVU632
.LBE199:
	.loc 1 21 42 is_stmt 1 discriminator 3 view .LVU633
	.loc 1 21 30 discriminator 1 view .LVU634
	testb	$1, %cl	#, niters.168
	je	.L108	#,
	.loc 1 21 30 is_stmt 0 discriminator 1 view .LVU635
	andl	$-2, %ecx	#, niters_vector_mult_vf.170
	addl	%ecx, %eax	# niters_vector_mult_vf.170,
.LVL179:
.L112:
.LBB200:
	.loc 1 22 9 is_stmt 1 view .LVU636
# src/kdtree.c:22:         FLOAT_TYPE dd = (p1[i] - p2[i]);
	.loc 1 22 20 is_stmt 0 view .LVU637
	vmovsd	(%rbx,%rax,8), %xmm1	# *_99, *_99
	vsubsd	(%rdx,%rax,8), %xmm1, %xmm1	# *_97, *_99, dd
.LVL180:
	.loc 1 23 9 is_stmt 1 view .LVU638
# src/kdtree.c:23:         d += dd*dd;
	.loc 1 23 16 is_stmt 0 view .LVU639
	vmulsd	%xmm1, %xmm1, %xmm1	# dd, dd, tmp177
.LVL181:
# src/kdtree.c:23:         d += dd*dd;
	.loc 1 23 11 view .LVU640
	vaddsd	%xmm1, %xmm0, %xmm0	# tmp177, d, d
.LVL182:
	.loc 1 23 11 view .LVU641
.LBE200:
	.loc 1 21 42 is_stmt 1 discriminator 3 view .LVU642
	.loc 1 21 30 discriminator 1 view .LVU643
.L108:
	.loc 1 21 30 is_stmt 0 discriminator 1 view .LVU644
.LBE202:
	.loc 1 25 2 is_stmt 1 view .LVU645
	.loc 1 25 2 is_stmt 0 view .LVU646
.LBE206:
.LBE210:
	.loc 1 226 5 is_stmt 1 view .LVU647
.LBB211:
.LBI211:
	.loc 1 212 19 view .LVU648
.LBB212:
	.loc 1 214 5 view .LVU649
# src/kdtree.c:214:     return p1[var] - p2[var];
	.loc 1 214 14 is_stmt 0 view .LVU650
	movslq	4(%r12), %rax	# kdtree_root_18->split_var, kdtree_root_18->split_var
.LBE212:
.LBE211:
# src/kdtree.c:227:     insertMaxHeap(H, current_distance, kdtree_root -> array_idx);
	.loc 1 227 5 view .LVU651
	movq	16(%r12), %rsi	# kdtree_root_18->array_idx, kdtree_root_18->array_idx
	movq	%r13, %rdi	# H,
.LBB214:
.LBB213:
# src/kdtree.c:214:     return p1[var] - p2[var];
	.loc 1 214 20 view .LVU652
	vmovsd	(%rbx,%rax,8), %xmm1	# *_30, *_30
	vsubsd	(%rdx,%rax,8), %xmm1, %xmm1	# *_32, *_30, _34
	vmovsd	%xmm1, -56(%rbp)	# _34, %sfp
.LVL183:
	.loc 1 214 20 view .LVU653
.LBE213:
.LBE214:
	.loc 1 227 5 is_stmt 1 view .LVU654
	call	insertMaxHeap@PLT	#
.LVL184:
	.loc 1 229 5 view .LVU655
	.loc 1 233 5 view .LVU656
	vmovsd	-56(%rbp), %xmm1	# %sfp, _34
	vxorpd	%xmm6, %xmm6, %xmm6	# tmp206
	vcomisd	%xmm6, %xmm1	# tmp206, _34
	jbe	.L140	#,
	.loc 1 240 13 view .LVU657
# src/kdtree.c:240:             if(kdtree_root -> rch) KNN_sub_tree_search(point, kdtree_root -> rch, H);
	.loc 1 240 28 is_stmt 0 view .LVU658
	movq	40(%r12), %rsi	# kdtree_root_18->rch, _7
# src/kdtree.c:240:             if(kdtree_root -> rch) KNN_sub_tree_search(point, kdtree_root -> rch, H);
	.loc 1 240 15 view .LVU659
	testq	%rsi, %rsi	# _7
	je	.L141	#,
# src/kdtree.c:240:             if(kdtree_root -> rch) KNN_sub_tree_search(point, kdtree_root -> rch, H);
	.loc 1 240 36 discriminator 1 view .LVU660
	movq	%r13, %rdx	# H,
	movq	%rbx, %rdi	# point,
	vmovsd	%xmm1, -56(%rbp)	# _34, %sfp
.LVL185:
	.loc 1 240 36 is_stmt 1 discriminator 1 view .LVU661
	call	KNN_sub_tree_search	#
.LVL186:
	.loc 1 246 5 view .LVU662
	.loc 1 247 5 view .LVU663
	.loc 1 250 5 view .LVU664
# src/kdtree.c:247:     int c   = max_d > (hp_distance * hp_distance);
	.loc 1 247 36 is_stmt 0 view .LVU665
	vmovsd	-56(%rbp), %xmm1	# %sfp, _34
# src/kdtree.c:246:     FLOAT_TYPE max_d = H -> data[0].value;
	.loc 1 246 16 view .LVU666
	movq	16(%r13), %rax	# H_22(D)->data, H_22(D)->data
# src/kdtree.c:247:     int c   = max_d > (hp_distance * hp_distance);
	.loc 1 247 36 view .LVU667
	vmulsd	%xmm1, %xmm1, %xmm1	# _34, _34, tmp193
# src/kdtree.c:250:     if(c || (H -> count) < (H -> N))
	.loc 1 250 7 view .LVU668
	vmovsd	(%rax), %xmm0	#* H_22(D)->data, _60->value
	vcomisd	%xmm1, %xmm0	# tmp193, _60->value
	ja	.L122	#,
.LVL187:
.L123:
# src/kdtree.c:250:     if(c || (H -> count) < (H -> N))
	.loc 1 250 10 discriminator 1 view .LVU669
	movq	0(%r13), %rax	# H_22(D)->N, tmp208
	cmpq	%rax, 8(%r13)	# tmp208, H_22(D)->count
	jnb	.L137	#,
.L122:
	.loc 1 260 17 is_stmt 1 view .LVU670
# src/kdtree.c:260:                 if(kdtree_root -> lch) KNN_sub_tree_search(point, kdtree_root -> lch, H);
	.loc 1 260 32 is_stmt 0 view .LVU671
	movq	32(%r12), %r12	# kdtree_root_18->lch, kdtree_root
# src/kdtree.c:260:                 if(kdtree_root -> lch) KNN_sub_tree_search(point, kdtree_root -> lch, H);
	.loc 1 260 19 view .LVU672
	testq	%r12, %r12	# kdtree_root
	jne	.L107	#,
.LVL188:
.L137:
# src/kdtree.c:268: }
	.loc 1 268 1 view .LVU673
	addq	$48, %rsp	#,
	popq	%rbx	#
.LVL189:
	.loc 1 268 1 view .LVU674
	popq	%r10	#
	.cfi_remember_state
	.cfi_def_cfa 10, 0
	popq	%r12	#
	popq	%r13	#
.LVL190:
	.loc 1 268 1 view .LVU675
	popq	%rbp	#
	leaq	-8(%r10), %rsp	#,
	.cfi_def_cfa 7, 8
	ret	
.LVL191:
	.p2align 4,,10
	.p2align 3
.L140:
	.cfi_restore_state
	.loc 1 236 13 is_stmt 1 view .LVU676
# src/kdtree.c:236:             if(kdtree_root -> lch) KNN_sub_tree_search(point, kdtree_root -> lch, H);
	.loc 1 236 28 is_stmt 0 view .LVU677
	movq	32(%r12), %rsi	# kdtree_root_18->lch, _6
# src/kdtree.c:236:             if(kdtree_root -> lch) KNN_sub_tree_search(point, kdtree_root -> lch, H);
	.loc 1 236 15 view .LVU678
	testq	%rsi, %rsi	# _6
	je	.L142	#,
# src/kdtree.c:236:             if(kdtree_root -> lch) KNN_sub_tree_search(point, kdtree_root -> lch, H);
	.loc 1 236 36 discriminator 1 view .LVU679
	movq	%r13, %rdx	# H,
	movq	%rbx, %rdi	# point,
	vmovsd	%xmm1, -56(%rbp)	# _34, %sfp
.LVL192:
	.loc 1 236 36 is_stmt 1 discriminator 1 view .LVU680
	call	KNN_sub_tree_search	#
.LVL193:
	.loc 1 246 5 view .LVU681
	.loc 1 247 5 view .LVU682
	.loc 1 250 5 view .LVU683
# src/kdtree.c:247:     int c   = max_d > (hp_distance * hp_distance);
	.loc 1 247 36 is_stmt 0 view .LVU684
	vmovsd	-56(%rbp), %xmm1	# %sfp, _34
# src/kdtree.c:246:     FLOAT_TYPE max_d = H -> data[0].value;
	.loc 1 246 16 view .LVU685
	movq	16(%r13), %rax	# H_22(D)->data, H_22(D)->data
# src/kdtree.c:247:     int c   = max_d > (hp_distance * hp_distance);
	.loc 1 247 36 view .LVU686
	vmulsd	%xmm1, %xmm1, %xmm1	# _34, _34, tmp187
# src/kdtree.c:250:     if(c || (H -> count) < (H -> N))
	.loc 1 250 7 view .LVU687
	vmovsd	(%rax), %xmm0	#* H_22(D)->data, _48->value
	vcomisd	%xmm1, %xmm0	# tmp187, _48->value
	ja	.L120	#,
.LVL194:
.L116:
# src/kdtree.c:250:     if(c || (H -> count) < (H -> N))
	.loc 1 250 10 discriminator 1 view .LVU688
	movq	0(%r13), %rax	# H_22(D)->N, tmp207
	cmpq	%rax, 8(%r13)	# tmp207, H_22(D)->count
	jnb	.L137	#,
.L120:
	.loc 1 256 17 is_stmt 1 view .LVU689
# src/kdtree.c:256:                 if(kdtree_root -> rch) KNN_sub_tree_search(point, kdtree_root -> rch, H);
	.loc 1 256 32 is_stmt 0 view .LVU690
	movq	40(%r12), %r12	# kdtree_root_18->rch, kdtree_root
# src/kdtree.c:256:                 if(kdtree_root -> rch) KNN_sub_tree_search(point, kdtree_root -> rch, H);
	.loc 1 256 19 view .LVU691
	testq	%r12, %r12	# kdtree_root
	je	.L137	#,
	.loc 1 225 5 is_stmt 1 view .LVU692
.LBB215:
.LBB207:
.LBB203:
# src/kdtree.c:21:     for(unsigned int i = 0; i<data_dims; ++i){
	.loc 1 21 30 is_stmt 0 discriminator 1 view .LVU693
	movl	data_dims(%rip), %ecx	# data_dims, data_dims.0_51
.LBE203:
.LBE207:
.LBE215:
# src/kdtree.c:225:     FLOAT_TYPE current_distance = euclidean_distance(point, kdtree_root -> data);
	.loc 1 225 35 view .LVU694
	movq	8(%r12), %rdx	# kdtree_root_18->data, _1
.LVL195:
.LBB216:
	.loc 1 19 12 is_stmt 1 view .LVU695
.LBB208:
	.loc 1 20 5 view .LVU696
	.loc 1 21 5 view .LVU697
.LBB204:
	.loc 1 21 9 view .LVU698
	.loc 1 21 30 discriminator 1 view .LVU699
	testl	%ecx, %ecx	# data_dims.0_51
	jne	.L143	#,
.LVL196:
.L124:
	.loc 1 21 30 is_stmt 0 discriminator 1 view .LVU700
.LBE204:
# src/kdtree.c:20:     FLOAT_TYPE d = 0;
	.loc 1 20 16 view .LVU701
	vxorpd	%xmm0, %xmm0, %xmm0	# d
	jmp	.L108	#
.LVL197:
	.p2align 4,,10
	.p2align 3
.L136:
	.loc 1 20 16 view .LVU702
	vzeroupper
.LVL198:
	.loc 1 20 16 view .LVU703
	jmp	.L108	#
.LVL199:
	.p2align 4,,10
	.p2align 3
.L141:
	.loc 1 20 16 view .LVU704
.LBE208:
.LBE216:
	.loc 1 246 5 is_stmt 1 view .LVU705
	.loc 1 247 5 view .LVU706
	.loc 1 250 5 view .LVU707
# src/kdtree.c:247:     int c   = max_d > (hp_distance * hp_distance);
	.loc 1 247 36 is_stmt 0 view .LVU708
	vmulsd	%xmm1, %xmm1, %xmm1	# _34, _34, tmp190
.LVL200:
# src/kdtree.c:246:     FLOAT_TYPE max_d = H -> data[0].value;
	.loc 1 246 16 view .LVU709
	movq	16(%r13), %rax	# H_22(D)->data, H_22(D)->data
# src/kdtree.c:250:     if(c || (H -> count) < (H -> N))
	.loc 1 250 7 view .LVU710
	vmovsd	(%rax), %xmm0	#* H_22(D)->data, _52->value
	vcomisd	%xmm1, %xmm0	# tmp190, _52->value
	ja	.L122	#,
	jmp	.L123	#
.LVL201:
	.p2align 4,,10
	.p2align 3
.L142:
	.loc 1 246 5 is_stmt 1 view .LVU711
	.loc 1 247 5 view .LVU712
	.loc 1 250 5 view .LVU713
# src/kdtree.c:247:     int c   = max_d > (hp_distance * hp_distance);
	.loc 1 247 36 is_stmt 0 view .LVU714
	vmulsd	%xmm1, %xmm1, %xmm1	# _34, _34, tmp184
.LVL202:
# src/kdtree.c:246:     FLOAT_TYPE max_d = H -> data[0].value;
	.loc 1 246 16 view .LVU715
	movq	16(%r13), %rax	# H_22(D)->data, H_22(D)->data
# src/kdtree.c:250:     if(c || (H -> count) < (H -> N))
	.loc 1 250 7 view .LVU716
	vmovsd	(%rax), %xmm0	#* H_22(D)->data, _56->value
	vcomisd	%xmm1, %xmm0	# tmp184, _56->value
	jbe	.L116	#,
	jmp	.L120	#
.LVL203:
.L125:
.LBB217:
.LBB209:
# src/kdtree.c:20:     FLOAT_TYPE d = 0;
	.loc 1 20 16 view .LVU717
	vxorpd	%xmm0, %xmm0, %xmm0	# d
.LBB205:
# src/kdtree.c:21:     for(unsigned int i = 0; i<data_dims; ++i){
	.loc 1 21 22 view .LVU718
	xorl	%eax, %eax	#
	jmp	.L109	#
.LBE205:
.LBE209:
.LBE217:
	.cfi_endproc
.LFE66:
	.size	KNN_sub_tree_search, .-KNN_sub_tree_search
	.p2align 4
	.globl	KNN
	.type	KNN, @function
KNN:
.LVL204:
.LFB67:
	.loc 1 274 1 is_stmt 1 view -0
	.cfi_startproc
	.loc 1 274 1 is_stmt 0 view .LVU720
	endbr64	
	pushq	%r13	#
	.cfi_def_cfa_offset 16
	.cfi_offset 13, -16
	movq	%rdx, %r13	# tmp98, kdtree_root
	pushq	%r12	#
	.cfi_def_cfa_offset 24
	.cfi_offset 12, -24
	movq	%rsi, %r12	# tmp97, point
# src/kdtree.c:276:     allocateHeap(&H,maxk);
	.loc 1 276 5 view .LVU721
	movslq	%ecx, %rsi	# tmp99, maxk
.LVL205:
# src/kdtree.c:274: {
	.loc 1 274 1 view .LVU722
	pushq	%rbp	#
	.cfi_def_cfa_offset 32
	.cfi_offset 6, -32
	pushq	%rbx	#
	.cfi_def_cfa_offset 40
	.cfi_offset 3, -40
	movq	%rdi, %rbx	# tmp96, .result_ptr
	subq	$40, %rsp	#,
	.cfi_def_cfa_offset 80
# src/kdtree.c:274: {
	.loc 1 274 1 view .LVU723
	movq	%fs:40, %rax	# MEM[(<address-space-1> long unsigned int *)40B], tmp100
	movq	%rax, 24(%rsp)	# tmp100, D.6650
	xorl	%eax, %eax	# tmp100
	.loc 1 275 5 is_stmt 1 view .LVU724
	.loc 1 276 5 view .LVU725
	movq	%rsp, %rbp	#, tmp89
	movq	%rbp, %rdi	# tmp89,
.LVL206:
	.loc 1 276 5 is_stmt 0 view .LVU726
	call	allocateHeap@PLT	#
.LVL207:
	.loc 1 277 5 is_stmt 1 view .LVU727
	movq	%rbp, %rdi	# tmp89,
	call	initHeap@PLT	#
.LVL208:
	.loc 1 278 5 view .LVU728
	movq	%rbp, %rdx	# tmp89,
	movq	%r13, %rsi	# kdtree_root,
	movq	%r12, %rdi	# point,
	call	KNN_sub_tree_search	#
.LVL209:
	.loc 1 279 5 view .LVU729
	movq	%rbp, %rdi	# tmp89,
	call	HeapSort@PLT	#
.LVL210:
	.loc 1 283 5 view .LVU730
# src/kdtree.c:283:     return H;
	.loc 1 283 12 is_stmt 0 view .LVU731
	movq	16(%rsp), %rax	# H, H
	vmovdqa	(%rsp), %xmm0	# H, tmp103
	movq	%rax, 16(%rbx)	# H, <retval>
	vmovdqu	%xmm0, (%rbx)	# tmp103, <retval>
# src/kdtree.c:284: }
	.loc 1 284 1 view .LVU732
	movq	24(%rsp), %rax	# D.6650, tmp101
	subq	%fs:40, %rax	# MEM[(<address-space-1> long unsigned int *)40B], tmp101
	jne	.L147	#,
	addq	$40, %rsp	#,
	.cfi_remember_state
	.cfi_def_cfa_offset 40
	movq	%rbx, %rax	# .result_ptr,
	popq	%rbx	#
	.cfi_def_cfa_offset 32
.LVL211:
	.loc 1 284 1 view .LVU733
	popq	%rbp	#
	.cfi_def_cfa_offset 24
	popq	%r12	#
	.cfi_def_cfa_offset 16
.LVL212:
	.loc 1 284 1 view .LVU734
	popq	%r13	#
	.cfi_def_cfa_offset 8
.LVL213:
	.loc 1 284 1 view .LVU735
	ret	
.LVL214:
.L147:
	.cfi_restore_state
	.loc 1 284 1 view .LVU736
	call	__stack_chk_fail@PLT	#
.LVL215:
	.cfi_endproc
.LFE67:
	.size	KNN, .-KNN
	.p2align 4
	.globl	build_tree
	.type	build_tree, @function
build_tree:
.LVL216:
.LFB68:
	.loc 1 287 1 is_stmt 1 view -0
	.cfi_startproc
	.loc 1 287 1 is_stmt 0 view .LVU738
	endbr64	
	.loc 1 296 5 is_stmt 1 view .LVU739
# src/kdtree.c:287: {
	.loc 1 287 1 is_stmt 0 view .LVU740
	pushq	%r13	#
	.cfi_def_cfa_offset 16
	.cfi_offset 13, -16
.LBB220:
.LBB221:
# src/kdtree.c:171:     if ((end - start) < 0) return 0;
	.loc 1 171 8 view .LVU741
	movl	%esi, %r13d	# tmp106, _4
.LBE221:
.LBE220:
# src/kdtree.c:287: {
	.loc 1 287 1 view .LVU742
	pushq	%r12	#
	.cfi_def_cfa_offset 24
	.cfi_offset 12, -24
	pushq	%rbp	#
	.cfi_def_cfa_offset 32
	.cfi_offset 6, -32
	pushq	%rbx	#
	.cfi_def_cfa_offset 40
	.cfi_offset 3, -40
	subq	$8, %rsp	#,
	.cfi_def_cfa_offset 48
# src/kdtree.c:296:    	data_dims = dimensions; 
	.loc 1 296 15 view .LVU743
	movl	%edx, data_dims(%rip)	# tmp107, data_dims
	.loc 1 298 5 is_stmt 1 view .LVU744
.LVL217:
.LBB227:
.LBI220:
	.loc 1 164 10 view .LVU745
.LBB222:
	.loc 1 166 5 view .LVU746
	.loc 1 167 5 view .LVU747
	.loc 1 170 5 view .LVU748
	.loc 1 171 5 view .LVU749
# src/kdtree.c:171:     if ((end - start) < 0) return 0;
	.loc 1 171 8 is_stmt 0 view .LVU750
	subl	$1, %r13d	#, _4
.LVL218:
	.loc 1 171 8 view .LVU751
	js	.L152	#,
	movq	%rdi, %rbp	# tmp105, kd_ptrs
	.loc 1 172 5 is_stmt 1 view .LVU752
# src/kdtree.c:172:     if (end  == start) {
	.loc 1 172 8 is_stmt 0 view .LVU753
	je	.L154	#,
	.loc 1 181 5 is_stmt 1 view .LVU754
# src/kdtree.c:181:     median_idx = medianOfNodes(t, start, end, split_var);
	.loc 1 181 18 is_stmt 0 view .LVU755
	xorl	%ecx, %ecx	#
	movl	%r13d, %edx	# _4,
.LVL219:
	.loc 1 181 18 view .LVU756
	xorl	%esi, %esi	#
.LVL220:
	.loc 1 181 18 view .LVU757
	call	medianOfNodes	#
.LVL221:
	.loc 1 181 18 view .LVU758
	movl	%eax, %r12d	# tmp108, median_idx
.LVL222:
	.loc 1 183 5 is_stmt 1 view .LVU759
# src/kdtree.c:183:     if(median_idx > -1){
	.loc 1 183 7 is_stmt 0 view .LVU760
	testl	%eax, %eax	# median_idx
	js	.L152	#,
	.loc 1 184 9 is_stmt 1 view .LVU761
# src/kdtree.c:184:         n = t[median_idx];
	.loc 1 184 14 is_stmt 0 view .LVU762
	cltq
.LVL223:
# src/kdtree.c:190: 				n->lch  = make_tree(t, start, median_idx - 1, n, level + 1);
	.loc 1 190 15 view .LVU763
	leal	-1(%r12), %edx	#, tmp101
	movq	%rbp, %rdi	# kd_ptrs,
	xorl	%esi, %esi	#
# src/kdtree.c:184:         n = t[median_idx];
	.loc 1 184 11 view .LVU764
	movq	0(%rbp,%rax,8), %rbx	# *_24, <retval>
.LVL224:
	.loc 1 190 5 is_stmt 1 view .LVU765
# src/kdtree.c:190: 				n->lch  = make_tree(t, start, median_idx - 1, n, level + 1);
	.loc 1 190 15 is_stmt 0 view .LVU766
	movq	%rbx, %rcx	# <retval>,
	call	make_tree.constprop.1	#
.LVL225:
# src/kdtree.c:192: 				n->rch = make_tree(t, median_idx + 1, end, n, level + 1);
	.loc 1 192 14 view .LVU767
	leal	1(%r12), %esi	#, tmp102
	movq	%rbx, %rcx	# <retval>,
	movl	%r13d, %edx	# _4,
# src/kdtree.c:190: 				n->lch  = make_tree(t, start, median_idx - 1, n, level + 1);
	.loc 1 190 13 discriminator 1 view .LVU768
	movq	%rax, 32(%rbx)	# tmp109, n_25->lch
	.loc 1 192 5 is_stmt 1 view .LVU769
# src/kdtree.c:192: 				n->rch = make_tree(t, median_idx + 1, end, n, level + 1);
	.loc 1 192 14 is_stmt 0 view .LVU770
	movq	%rbp, %rdi	# kd_ptrs,
	call	make_tree.constprop.1	#
.LVL226:
# src/kdtree.c:196:         n->parent = parent;
	.loc 1 196 19 view .LVU771
	movq	$0, 24(%rbx)	#, n_25->parent
# src/kdtree.c:192: 				n->rch = make_tree(t, median_idx + 1, end, n, level + 1);
	.loc 1 192 12 discriminator 1 view .LVU772
	movq	%rax, 40(%rbx)	# tmp110, n_25->rch
	.loc 1 195 9 is_stmt 1 view .LVU773
	.loc 1 196 9 view .LVU774
	.loc 1 197 9 view .LVU775
.LBE222:
.LBE227:
# src/kdtree.c:302: }
	.loc 1 302 1 is_stmt 0 view .LVU776
	movq	%rbx, %rax	# <retval>,
.LBB228:
.LBB223:
# src/kdtree.c:197:         n->level = level;
	.loc 1 197 18 view .LVU777
	movq	$0, (%rbx)	#, MEM <vector(2) int> [(int *)n_25]
.LVL227:
	.loc 1 197 18 view .LVU778
.LBE223:
.LBE228:
	.loc 1 300 5 is_stmt 1 view .LVU779
# src/kdtree.c:302: }
	.loc 1 302 1 is_stmt 0 view .LVU780
	addq	$8, %rsp	#,
	.cfi_remember_state
	.cfi_def_cfa_offset 40
	popq	%rbx	#
	.cfi_def_cfa_offset 32
.LVL228:
	.loc 1 302 1 view .LVU781
	popq	%rbp	#
	.cfi_def_cfa_offset 24
.LVL229:
	.loc 1 302 1 view .LVU782
	popq	%r12	#
	.cfi_def_cfa_offset 16
	popq	%r13	#
	.cfi_def_cfa_offset 8
	ret	
.LVL230:
	.p2align 4,,10
	.p2align 3
.L152:
	.cfi_restore_state
	.loc 1 302 1 view .LVU783
	addq	$8, %rsp	#,
	.cfi_remember_state
	.cfi_def_cfa_offset 40
.LBB229:
.LBB224:
# src/kdtree.c:171:     if ((end - start) < 0) return 0;
	.loc 1 171 35 discriminator 1 view .LVU784
	xorl	%ebx, %ebx	# <retval>
.LBE224:
.LBE229:
# src/kdtree.c:302: }
	.loc 1 302 1 view .LVU785
	movq	%rbx, %rax	# <retval>,
	popq	%rbx	#
	.cfi_def_cfa_offset 32
	popq	%rbp	#
	.cfi_def_cfa_offset 24
	popq	%r12	#
	.cfi_def_cfa_offset 16
	popq	%r13	#
	.cfi_def_cfa_offset 8
.LVL231:
	.loc 1 302 1 view .LVU786
	ret	
.LVL232:
	.p2align 4,,10
	.p2align 3
.L154:
	.cfi_restore_state
.LBB230:
.LBB225:
	.loc 1 173 9 is_stmt 1 view .LVU787
# src/kdtree.c:173:         n = t[start];
	.loc 1 173 11 is_stmt 0 view .LVU788
	movq	(%rdi), %rbx	# *kd_ptrs_9(D), <retval>
.LVL233:
	.loc 1 174 9 is_stmt 1 view .LVU789
	.loc 1 175 9 view .LVU790
	.loc 1 176 9 view .LVU791
# src/kdtree.c:175:         n->parent = parent;
	.loc 1 175 19 is_stmt 0 view .LVU792
	vpxor	%xmm0, %xmm0, %xmm0	# tmp99
# src/kdtree.c:176:         n->level = level;
	.loc 1 176 18 view .LVU793
	movq	$0, (%rbx)	#, MEM <vector(2) int> [(int *)n_20]
	.loc 1 177 9 is_stmt 1 view .LVU794
.LBE225:
.LBE230:
# src/kdtree.c:302: }
	.loc 1 302 1 is_stmt 0 view .LVU795
	movq	%rbx, %rax	# <retval>,
.LBB231:
.LBB226:
# src/kdtree.c:178:         n -> rch = NULL;
	.loc 1 178 18 view .LVU796
	movq	$0, 40(%rbx)	#, n_20->rch
# src/kdtree.c:175:         n->parent = parent;
	.loc 1 175 19 view .LVU797
	vmovdqu	%xmm0, 24(%rbx)	# tmp99, MEM <vector(2) long unsigned int> [(struct kd_node * *)n_20 + 24B]
	.loc 1 178 9 is_stmt 1 view .LVU798
	.loc 1 179 9 view .LVU799
.LBE226:
.LBE231:
# src/kdtree.c:302: }
	.loc 1 302 1 is_stmt 0 view .LVU800
	addq	$8, %rsp	#,
	.cfi_def_cfa_offset 40
	popq	%rbx	#
	.cfi_def_cfa_offset 32
.LVL234:
	.loc 1 302 1 view .LVU801
	popq	%rbp	#
	.cfi_def_cfa_offset 24
	popq	%r12	#
	.cfi_def_cfa_offset 16
	popq	%r13	#
	.cfi_def_cfa_offset 8
.LVL235:
	.loc 1 302 1 view .LVU802
	ret	
	.cfi_endproc
.LFE68:
	.size	build_tree, .-build_tree
	.globl	c
	.bss
	.align 4
	.type	c, @object
	.size	c, 4
c:
	.zero	4
	.section	.rodata.cst8,"aM",@progbits,8
	.align 8
.LC1:
	.long	-1
	.long	-1
	.section	.rodata.cst32,"aM",@progbits,32
	.align 32
.LC3:
	.quad	0
	.quad	1
	.quad	2
	.quad	3
	.section	.rodata.cst8
	.align 8
.LC5:
	.quad	4
	.text
.Letext0:
	.file 4 "/usr/lib/gcc/x86_64-linux-gnu/13/include/stddef.h"
	.file 5 "/usr/include/x86_64-linux-gnu/bits/types.h"
	.file 6 "/usr/include/x86_64-linux-gnu/bits/stdint-uintn.h"
	.file 7 "src/../include/heap.h"
	.file 8 "src/../include/kdtree.h"
	.file 9 "/usr/include/x86_64-linux-gnu/bits/stdio2-decl.h"
	.file 10 "<built-in>"
	.section	.debug_info,"",@progbits
.Ldebug_info0:
	.long	0x17e9
	.value	0x5
	.byte	0x1
	.byte	0x8
	.long	.Ldebug_abbrev0
	.uleb128 0x2c
	.long	.LASF73
	.byte	0x1d
	.long	.LASF0
	.long	.LASF1
	.quad	.Ltext0
	.quad	.Letext0-.Ltext0
	.long	.Ldebug_line0
	.uleb128 0x7
	.byte	0x8
	.byte	0x4
	.long	.LASF2
	.uleb128 0x11
	.long	.LASF10
	.byte	0x4
	.byte	0xd6
	.byte	0x17
	.long	0x41
	.uleb128 0x7
	.byte	0x8
	.byte	0x7
	.long	.LASF3
	.uleb128 0x7
	.byte	0x4
	.byte	0x7
	.long	.LASF4
	.uleb128 0x2d
	.byte	0x8
	.uleb128 0x22
	.long	0x4f
	.uleb128 0x7
	.byte	0x1
	.byte	0x8
	.long	.LASF5
	.uleb128 0x7
	.byte	0x2
	.byte	0x7
	.long	.LASF6
	.uleb128 0x7
	.byte	0x1
	.byte	0x6
	.long	.LASF7
	.uleb128 0x7
	.byte	0x2
	.byte	0x5
	.long	.LASF8
	.uleb128 0x2e
	.byte	0x4
	.byte	0x5
	.string	"int"
	.uleb128 0x7
	.byte	0x8
	.byte	0x5
	.long	.LASF9
	.uleb128 0x11
	.long	.LASF11
	.byte	0x5
	.byte	0x2d
	.byte	0x1b
	.long	0x41
	.uleb128 0x7
	.byte	0x1
	.byte	0x6
	.long	.LASF12
	.uleb128 0x2f
	.long	0x8c
	.uleb128 0xb
	.long	0x93
	.uleb128 0x22
	.long	0x98
	.uleb128 0x7
	.byte	0x8
	.byte	0x5
	.long	.LASF13
	.uleb128 0x7
	.byte	0x8
	.byte	0x7
	.long	.LASF14
	.uleb128 0xb
	.long	0xba
	.uleb128 0x22
	.long	0xb0
	.uleb128 0x30
	.uleb128 0x7
	.byte	0x4
	.byte	0x4
	.long	.LASF15
	.uleb128 0x11
	.long	.LASF16
	.byte	0x6
	.byte	0x1b
	.byte	0x14
	.long	0x80
	.uleb128 0x23
	.long	.LASF19
	.byte	0x10
	.byte	0x7
	.byte	0x22
	.long	0xf5
	.uleb128 0xc
	.long	.LASF17
	.byte	0x7
	.byte	0x24
	.byte	0xf
	.long	0x2e
	.byte	0
	.uleb128 0xc
	.long	.LASF18
	.byte	0x7
	.byte	0x25
	.byte	0xa
	.long	0xc2
	.byte	0x8
	.byte	0
	.uleb128 0x23
	.long	.LASF20
	.byte	0x18
	.byte	0x7
	.byte	0x28
	.long	0x127
	.uleb128 0x24
	.string	"N"
	.byte	0x7
	.byte	0x2a
	.byte	0xa
	.long	0xc2
	.byte	0
	.uleb128 0xc
	.long	.LASF21
	.byte	0x7
	.byte	0x2b
	.byte	0xa
	.long	0xc2
	.byte	0x8
	.uleb128 0xc
	.long	.LASF22
	.byte	0x7
	.byte	0x2c
	.byte	0x16
	.long	0x127
	.byte	0x10
	.byte	0
	.uleb128 0xb
	.long	0xce
	.uleb128 0xb
	.long	0x2e
	.uleb128 0x11
	.long	.LASF20
	.byte	0x7
	.byte	0x37
	.byte	0x15
	.long	0xf5
	.uleb128 0x11
	.long	.LASF19
	.byte	0x7
	.byte	0x38
	.byte	0x1a
	.long	0xce
	.uleb128 0x23
	.long	.LASF23
	.byte	0x30
	.byte	0x8
	.byte	0x1c
	.long	0x1b1
	.uleb128 0xc
	.long	.LASF24
	.byte	0x8
	.byte	0x1e
	.byte	0x8
	.long	0x72
	.byte	0
	.uleb128 0xc
	.long	.LASF25
	.byte	0x8
	.byte	0x1f
	.byte	0x8
	.long	0x72
	.byte	0x4
	.uleb128 0xc
	.long	.LASF22
	.byte	0x8
	.byte	0x20
	.byte	0x11
	.long	0x12c
	.byte	0x8
	.uleb128 0xc
	.long	.LASF18
	.byte	0x8
	.byte	0x21
	.byte	0xa
	.long	0xc2
	.byte	0x10
	.uleb128 0xc
	.long	.LASF26
	.byte	0x8
	.byte	0x22
	.byte	0x14
	.long	0x1b1
	.byte	0x18
	.uleb128 0x24
	.string	"lch"
	.byte	0x8
	.byte	0x23
	.byte	0x14
	.long	0x1b1
	.byte	0x20
	.uleb128 0x24
	.string	"rch"
	.byte	0x8
	.byte	0x24
	.byte	0x14
	.long	0x1b1
	.byte	0x28
	.byte	0
	.uleb128 0xb
	.long	0x149
	.uleb128 0x11
	.long	.LASF23
	.byte	0x8
	.byte	0x27
	.byte	0x18
	.long	0x149
	.uleb128 0x31
	.long	.LASF27
	.byte	0x1
	.byte	0x8
	.byte	0x15
	.long	0x48
	.uleb128 0x32
	.string	"c"
	.byte	0x1
	.byte	0xa
	.byte	0x5
	.long	0x72
	.uleb128 0x9
	.byte	0x3
	.quad	c
	.uleb128 0x1b
	.long	.LASF28
	.byte	0x50
	.long	0x1f2
	.uleb128 0xd
	.long	0x1f2
	.byte	0
	.uleb128 0xb
	.long	0x131
	.uleb128 0x1b
	.long	.LASF29
	.byte	0x40
	.long	0x207
	.uleb128 0xd
	.long	0x1f2
	.byte	0
	.uleb128 0x1b
	.long	.LASF30
	.byte	0x3c
	.long	0x21c
	.uleb128 0xd
	.long	0x1f2
	.uleb128 0xd
	.long	0xc2
	.byte	0
	.uleb128 0x1b
	.long	.LASF31
	.byte	0x4e
	.long	0x236
	.uleb128 0xd
	.long	0x1f2
	.uleb128 0xd
	.long	0x2e
	.uleb128 0xd
	.long	0xc2
	.byte	0
	.uleb128 0x33
	.long	.LASF43
	.byte	0x9
	.byte	0x33
	.byte	0xc
	.long	0x72
	.long	0x252
	.uleb128 0xd
	.long	0x72
	.uleb128 0xd
	.long	0x98
	.uleb128 0x28
	.byte	0
	.uleb128 0x34
	.long	.LASF35
	.byte	0x1
	.value	0x11e
	.byte	0xb
	.long	0x3d0
	.quad	.LFB68
	.quad	.LFE68-.LFB68
	.uleb128 0x1
	.byte	0x9c
	.long	0x3d0
	.uleb128 0x17
	.long	.LASF32
	.value	0x11e
	.byte	0x20
	.long	0x3d5
	.long	.LLST145
	.long	.LVUS145
	.uleb128 0x35
	.string	"n"
	.byte	0x1
	.value	0x11e
	.byte	0x30
	.long	0x35
	.long	.LLST146
	.long	.LVUS146
	.uleb128 0x17
	.long	.LASF33
	.value	0x11e
	.byte	0x3a
	.long	0x35
	.long	.LLST147
	.long	.LVUS147
	.uleb128 0x36
	.long	.LASF34
	.byte	0x1
	.value	0x12a
	.byte	0xe
	.long	0x3d0
	.long	.LLST148
	.long	.LVUS148
	.uleb128 0x37
	.long	0x6c6
	.quad	.LBI220
	.byte	.LVU745
	.long	.LLRL149
	.byte	0x1
	.value	0x12a
	.byte	0x15
	.uleb128 0x2
	.long	0x6df
	.long	.LLST150
	.long	.LVUS150
	.uleb128 0x2
	.long	0x6f6
	.long	.LLST150
	.long	.LVUS150
	.uleb128 0x2
	.long	0x702
	.long	.LLST150
	.long	.LVUS150
	.uleb128 0x2
	.long	0x6eb
	.long	.LLST153
	.long	.LVUS153
	.uleb128 0x2
	.long	0x6d6
	.long	.LLST154
	.long	.LVUS154
	.uleb128 0xe
	.long	.LLRL149
	.uleb128 0x3
	.long	0x70e
	.long	.LLST155
	.long	.LVUS155
	.uleb128 0x3
	.long	0x717
	.long	.LLST156
	.long	.LVUS156
	.uleb128 0x3
	.long	0x722
	.long	.LLST157
	.long	.LVUS157
	.uleb128 0x4
	.quad	.LVL221
	.long	0x72e
	.long	0x370
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x55
	.uleb128 0x2
	.byte	0x76
	.sleb128 0
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x54
	.uleb128 0x1
	.byte	0x30
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x51
	.uleb128 0x2
	.byte	0x7d
	.sleb128 0
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x52
	.uleb128 0x1
	.byte	0x30
	.byte	0
	.uleb128 0x4
	.quad	.LVL225
	.long	0x16c3
	.long	0x3a0
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x55
	.uleb128 0x2
	.byte	0x76
	.sleb128 0
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x54
	.uleb128 0x1
	.byte	0x30
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x51
	.uleb128 0x2
	.byte	0x7c
	.sleb128 -1
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x52
	.uleb128 0x2
	.byte	0x73
	.sleb128 0
	.uleb128 0x1c
	.long	0x702
	.uleb128 0x1
	.byte	0x31
	.byte	0
	.uleb128 0x6
	.quad	.LVL226
	.long	0x16c3
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x55
	.uleb128 0x2
	.byte	0x76
	.sleb128 0
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x54
	.uleb128 0x2
	.byte	0x7c
	.sleb128 1
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x51
	.uleb128 0x2
	.byte	0x7d
	.sleb128 0
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x52
	.uleb128 0x2
	.byte	0x73
	.sleb128 0
	.uleb128 0x1c
	.long	0x702
	.uleb128 0x1
	.byte	0x31
	.byte	0
	.byte	0
	.byte	0
	.byte	0
	.uleb128 0xb
	.long	0x1b6
	.uleb128 0xb
	.long	0x3d0
	.uleb128 0x38
	.string	"KNN"
	.byte	0x1
	.value	0x111
	.byte	0x6
	.long	0x131
	.quad	.LFB67
	.quad	.LFE67-.LFB67
	.uleb128 0x1
	.byte	0x9c
	.long	0x4cf
	.uleb128 0x17
	.long	.LASF36
	.value	0x111
	.byte	0x16
	.long	0x12c
	.long	.LLST142
	.long	.LVUS142
	.uleb128 0x17
	.long	.LASF37
	.value	0x111
	.byte	0x26
	.long	0x3d0
	.long	.LLST143
	.long	.LVUS143
	.uleb128 0x17
	.long	.LASF38
	.value	0x111
	.byte	0x37
	.long	0x72
	.long	.LLST144
	.long	.LVUS144
	.uleb128 0x39
	.string	"H"
	.byte	0x1
	.value	0x113
	.byte	0xa
	.long	0x131
	.uleb128 0x3
	.byte	0x91
	.sleb128 -80
	.uleb128 0x4
	.quad	.LVL207
	.long	0x207
	.long	0x46d
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x55
	.uleb128 0x2
	.byte	0x76
	.sleb128 0
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x54
	.uleb128 0x9
	.byte	0xa3
	.uleb128 0x1
	.byte	0x52
	.byte	0x8
	.byte	0x20
	.byte	0x24
	.byte	0x8
	.byte	0x20
	.byte	0x26
	.byte	0
	.uleb128 0x4
	.quad	.LVL208
	.long	0x1f7
	.long	0x485
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x55
	.uleb128 0x2
	.byte	0x76
	.sleb128 0
	.byte	0
	.uleb128 0x4
	.quad	.LVL209
	.long	0x4cf
	.long	0x4a9
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x55
	.uleb128 0x2
	.byte	0x7c
	.sleb128 0
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x54
	.uleb128 0x2
	.byte	0x7d
	.sleb128 0
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x51
	.uleb128 0x2
	.byte	0x76
	.sleb128 0
	.byte	0
	.uleb128 0x4
	.quad	.LVL210
	.long	0x1e2
	.long	0x4c1
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x55
	.uleb128 0x2
	.byte	0x76
	.sleb128 0
	.byte	0
	.uleb128 0x3a
	.quad	.LVL215
	.long	0x17d2
	.byte	0
	.uleb128 0x12
	.long	.LASF55
	.byte	0xde
	.quad	.LFB66
	.quad	.LFE66-.LFB66
	.uleb128 0x1
	.byte	0x9c
	.long	0x696
	.uleb128 0x1d
	.long	.LASF36
	.byte	0xde
	.byte	0x26
	.long	0x12c
	.long	.LLST123
	.long	.LVUS123
	.uleb128 0x1d
	.long	.LASF37
	.byte	0xde
	.byte	0x36
	.long	0x3d0
	.long	.LLST124
	.long	.LVUS124
	.uleb128 0x25
	.string	"H"
	.byte	0xde
	.byte	0x4a
	.long	0x1f2
	.long	.LLST125
	.long	.LVUS125
	.uleb128 0x1e
	.long	.LASF39
	.byte	0xe1
	.byte	0x10
	.long	0x2e
	.long	.LLST126
	.long	.LVUS126
	.uleb128 0x1e
	.long	.LASF40
	.byte	0xe2
	.byte	0x10
	.long	0x2e
	.long	.LLST127
	.long	.LVUS127
	.uleb128 0x1e
	.long	.LASF41
	.byte	0xe5
	.byte	0x9
	.long	0x72
	.long	.LLST128
	.long	.LVUS128
	.uleb128 0x1e
	.long	.LASF42
	.byte	0xf6
	.byte	0x10
	.long	0x2e
	.long	.LLST129
	.long	.LVUS129
	.uleb128 0x1f
	.string	"c"
	.byte	0xf7
	.byte	0x9
	.long	0x72
	.long	.LLST130
	.long	.LVUS130
	.uleb128 0x13
	.long	0xdca
	.quad	.LBI195
	.byte	.LVU609
	.long	.LLRL131
	.byte	0xe1
	.byte	0x23
	.long	0x5eb
	.uleb128 0x2
	.long	0xde4
	.long	.LLST132
	.long	.LVUS132
	.uleb128 0x18
	.long	0xdda
	.uleb128 0xe
	.long	.LLRL131
	.uleb128 0x3
	.long	0xdee
	.long	.LLST133
	.long	.LVUS133
	.uleb128 0x26
	.long	0xdf7
	.long	.LLRL134
	.uleb128 0x3
	.long	0xdf8
	.long	.LLST135
	.long	.LVUS135
	.uleb128 0x26
	.long	0xe01
	.long	.LLRL136
	.uleb128 0x3
	.long	0xe02
	.long	.LLST137
	.long	.LVUS137
	.byte	0
	.byte	0
	.byte	0
	.byte	0
	.uleb128 0x13
	.long	0x696
	.quad	.LBI211
	.byte	.LVU648
	.long	.LLRL138
	.byte	0xe2
	.byte	0x1e
	.long	0x62b
	.uleb128 0x2
	.long	0x6ba
	.long	.LLST139
	.long	.LVUS139
	.uleb128 0x2
	.long	0x6b0
	.long	.LLST140
	.long	.LVUS140
	.uleb128 0x2
	.long	0x6a6
	.long	.LLST141
	.long	.LVUS141
	.byte	0
	.uleb128 0x14
	.quad	.LVL174
	.uleb128 0x1
	.byte	0x30
	.uleb128 0x4
	.quad	.LVL184
	.long	0x21c
	.long	0x64e
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x55
	.uleb128 0x2
	.byte	0x7d
	.sleb128 0
	.byte	0
	.uleb128 0x4
	.quad	.LVL186
	.long	0x4cf
	.long	0x66c
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x55
	.uleb128 0x2
	.byte	0x73
	.sleb128 0
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x51
	.uleb128 0x2
	.byte	0x7d
	.sleb128 0
	.byte	0
	.uleb128 0x4
	.quad	.LVL193
	.long	0x4cf
	.long	0x68a
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x55
	.uleb128 0x2
	.byte	0x73
	.sleb128 0
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x51
	.uleb128 0x2
	.byte	0x7d
	.sleb128 0
	.byte	0
	.uleb128 0x14
	.quad	.LVL198
	.uleb128 0x1
	.byte	0x30
	.byte	0
	.uleb128 0x15
	.long	.LASF44
	.byte	0xd4
	.byte	0x13
	.long	0x2e
	.byte	0x3
	.long	0x6c6
	.uleb128 0x5
	.string	"p1"
	.byte	0xd4
	.byte	0x30
	.long	0x12c
	.uleb128 0x5
	.string	"p2"
	.byte	0xd4
	.byte	0x40
	.long	0x12c
	.uleb128 0x5
	.string	"var"
	.byte	0xd4
	.byte	0x48
	.long	0x72
	.byte	0
	.uleb128 0x15
	.long	.LASF45
	.byte	0xa4
	.byte	0xa
	.long	0x3d0
	.byte	0x1
	.long	0x72e
	.uleb128 0x5
	.string	"t"
	.byte	0xa4
	.byte	0x1e
	.long	0x3d5
	.uleb128 0x8
	.long	.LASF46
	.byte	0x1
	.byte	0xa4
	.byte	0x25
	.long	0x72
	.uleb128 0x5
	.string	"end"
	.byte	0xa4
	.byte	0x30
	.long	0x72
	.uleb128 0x8
	.long	.LASF26
	.byte	0x1
	.byte	0xa4
	.byte	0x3e
	.long	0x3d0
	.uleb128 0x8
	.long	.LASF24
	.byte	0x1
	.byte	0xa4
	.byte	0x4a
	.long	0x72
	.uleb128 0xa
	.string	"n"
	.byte	0xa6
	.byte	0xe
	.long	0x3d0
	.uleb128 0x20
	.long	.LASF25
	.byte	0xa7
	.byte	0x9
	.long	0x72
	.uleb128 0x20
	.long	.LASF47
	.byte	0xaa
	.byte	0x9
	.long	0x72
	.byte	0
	.uleb128 0x15
	.long	.LASF48
	.byte	0x74
	.byte	0x5
	.long	0x72
	.byte	0x1
	.long	0x782
	.uleb128 0x5
	.string	"a"
	.byte	0x74
	.byte	0x1d
	.long	0x3d5
	.uleb128 0x8
	.long	.LASF49
	.byte	0x1
	.byte	0x74
	.byte	0x24
	.long	0x72
	.uleb128 0x8
	.long	.LASF50
	.byte	0x1
	.byte	0x74
	.byte	0x2e
	.long	0x72
	.uleb128 0x8
	.long	.LASF25
	.byte	0x1
	.byte	0x74
	.byte	0x39
	.long	0x72
	.uleb128 0xa
	.string	"k"
	.byte	0x77
	.byte	0x9
	.long	0x72
	.uleb128 0x21
	.uleb128 0x20
	.long	.LASF51
	.byte	0x8f
	.byte	0xd
	.long	0x72
	.byte	0
	.byte	0
	.uleb128 0x15
	.long	.LASF52
	.byte	0x64
	.byte	0x5
	.long	0x72
	.byte	0x1
	.long	0x7e0
	.uleb128 0x5
	.string	"arr"
	.byte	0x64
	.byte	0x19
	.long	0x3d5
	.uleb128 0x5
	.string	"low"
	.byte	0x64
	.byte	0x22
	.long	0x72
	.uleb128 0x8
	.long	.LASF53
	.byte	0x1
	.byte	0x64
	.byte	0x2b
	.long	0x72
	.uleb128 0x8
	.long	.LASF25
	.byte	0x1
	.byte	0x64
	.byte	0x35
	.long	0x72
	.uleb128 0x20
	.long	.LASF54
	.byte	0x66
	.byte	0xe
	.long	0x3d0
	.uleb128 0xa
	.string	"i"
	.byte	0x68
	.byte	0x9
	.long	0x72
	.uleb128 0x21
	.uleb128 0xa
	.string	"j"
	.byte	0x69
	.byte	0xe
	.long	0x72
	.byte	0
	.byte	0
	.uleb128 0x12
	.long	.LASF56
	.byte	0x53
	.quad	.LFB60
	.quad	.LFE60-.LFB60
	.uleb128 0x1
	.byte	0x9c
	.long	0xb6d
	.uleb128 0x1d
	.long	.LASF57
	.byte	0x53
	.byte	0x1b
	.long	0x3d0
	.long	.LLST30
	.long	.LVUS30
	.uleb128 0x3b
	.quad	.LBB95
	.quad	.LBE95-.LBB95
	.long	0x871
	.uleb128 0x1f
	.string	"i"
	.byte	0x58
	.byte	0x16
	.long	0x48
	.long	.LLST35
	.long	.LVUS35
	.uleb128 0x19
	.long	0xf5f
	.quad	.LBI96
	.byte	.LVU195
	.long	.LLRL36
	.byte	0x58
	.byte	0x2d
	.uleb128 0x2
	.long	0xf6e
	.long	.LLST37
	.long	.LVUS37
	.uleb128 0x6
	.quad	.LVL56
	.long	0x236
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x55
	.uleb128 0x1
	.byte	0x31
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x54
	.uleb128 0x2
	.byte	0x7c
	.sleb128 0
	.byte	0
	.byte	0
	.byte	0
	.uleb128 0x13
	.long	0xf5f
	.quad	.LBI83
	.byte	.LVU171
	.long	.LLRL31
	.byte	0x55
	.byte	0x5
	.long	0x8bd
	.uleb128 0x2
	.long	0xf6e
	.long	.LLST32
	.long	.LVUS32
	.uleb128 0x6
	.quad	.LVL51
	.long	0x236
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x55
	.uleb128 0x1
	.byte	0x31
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x54
	.uleb128 0x9
	.byte	0x3
	.quad	.LC6
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x51
	.uleb128 0x2
	.byte	0x76
	.sleb128 0
	.byte	0
	.byte	0
	.uleb128 0x9
	.long	0xf5f
	.quad	.LBI91
	.byte	.LVU182
	.quad	.LBB91
	.quad	.LBE91-.LBB91
	.byte	0x56
	.long	0x90e
	.uleb128 0x2
	.long	0xf6e
	.long	.LLST33
	.long	.LVUS33
	.uleb128 0x6
	.quad	.LVL52
	.long	0x236
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x55
	.uleb128 0x1
	.byte	0x31
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x54
	.uleb128 0x9
	.byte	0x3
	.quad	.LC7
	.byte	0
	.byte	0
	.uleb128 0x9
	.long	0xf5f
	.quad	.LBI93
	.byte	.LVU187
	.quad	.LBB93
	.quad	.LBE93-.LBB93
	.byte	0x57
	.long	0x95f
	.uleb128 0x2
	.long	0xf6e
	.long	.LLST34
	.long	.LVUS34
	.uleb128 0x6
	.quad	.LVL53
	.long	0x236
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x55
	.uleb128 0x1
	.byte	0x31
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x54
	.uleb128 0x9
	.byte	0x3
	.quad	.LC8
	.byte	0
	.byte	0
	.uleb128 0x9
	.long	0xf5f
	.quad	.LBI102
	.byte	.LVU206
	.quad	.LBB102
	.quad	.LBE102-.LBB102
	.byte	0x59
	.long	0x9a3
	.uleb128 0x2
	.long	0xf6e
	.long	.LLST38
	.long	.LVUS38
	.uleb128 0x6
	.quad	.LVL58
	.long	0x17e1
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x55
	.uleb128 0x1
	.byte	0x3a
	.byte	0
	.byte	0
	.uleb128 0x9
	.long	0xf5f
	.quad	.LBI104
	.byte	.LVU211
	.quad	.LBB104
	.quad	.LBE104-.LBB104
	.byte	0x5a
	.long	0x9f4
	.uleb128 0x2
	.long	0xf6e
	.long	.LLST39
	.long	.LVUS39
	.uleb128 0x6
	.quad	.LVL59
	.long	0x236
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x55
	.uleb128 0x1
	.byte	0x31
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x54
	.uleb128 0x9
	.byte	0x3
	.quad	.LC10
	.byte	0
	.byte	0
	.uleb128 0x9
	.long	0xf5f
	.quad	.LBI106
	.byte	.LVU216
	.quad	.LBB106
	.quad	.LBE106-.LBB106
	.byte	0x5b
	.long	0xa45
	.uleb128 0x2
	.long	0xf6e
	.long	.LLST40
	.long	.LVUS40
	.uleb128 0x6
	.quad	.LVL60
	.long	0x236
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x55
	.uleb128 0x1
	.byte	0x31
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x54
	.uleb128 0x9
	.byte	0x3
	.quad	.LC11
	.byte	0
	.byte	0
	.uleb128 0x9
	.long	0xf5f
	.quad	.LBI108
	.byte	.LVU221
	.quad	.LBB108
	.quad	.LBE108-.LBB108
	.byte	0x5c
	.long	0xa96
	.uleb128 0x2
	.long	0xf6e
	.long	.LLST41
	.long	.LVUS41
	.uleb128 0x6
	.quad	.LVL61
	.long	0x236
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x55
	.uleb128 0x1
	.byte	0x31
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x54
	.uleb128 0x9
	.byte	0x3
	.quad	.LC12
	.byte	0
	.byte	0
	.uleb128 0x9
	.long	0xf5f
	.quad	.LBI110
	.byte	.LVU226
	.quad	.LBB110
	.quad	.LBE110-.LBB110
	.byte	0x5d
	.long	0xae7
	.uleb128 0x2
	.long	0xf6e
	.long	.LLST42
	.long	.LVUS42
	.uleb128 0x6
	.quad	.LVL62
	.long	0x236
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x55
	.uleb128 0x1
	.byte	0x31
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x54
	.uleb128 0x9
	.byte	0x3
	.quad	.LC13
	.byte	0
	.byte	0
	.uleb128 0x9
	.long	0xf5f
	.quad	.LBI112
	.byte	.LVU231
	.quad	.LBB112
	.quad	.LBE112-.LBB112
	.byte	0x5e
	.long	0xb38
	.uleb128 0x2
	.long	0xf6e
	.long	.LLST43
	.long	.LVUS43
	.uleb128 0x6
	.quad	.LVL63
	.long	0x236
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x55
	.uleb128 0x1
	.byte	0x31
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x54
	.uleb128 0x9
	.byte	0x3
	.quad	.LC14
	.byte	0
	.byte	0
	.uleb128 0x19
	.long	0xf5f
	.quad	.LBI114
	.byte	.LVU236
	.long	.LLRL44
	.byte	0x5f
	.byte	0x5
	.uleb128 0xf
	.long	0xf6e
	.uleb128 0x6
	.byte	0xa0
	.long	.Ldebug_info0+6107
	.sleb128 0
	.uleb128 0x3c
	.quad	.LVL65
	.long	0x17e1
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x55
	.uleb128 0x1
	.byte	0x3a
	.byte	0
	.byte	0
	.byte	0
	.uleb128 0x15
	.long	.LASF58
	.byte	0x4d
	.byte	0x5
	.long	0x72
	.byte	0x1
	.long	0xba6
	.uleb128 0x5
	.string	"a"
	.byte	0x4d
	.byte	0x19
	.long	0x3d0
	.uleb128 0x5
	.string	"b"
	.byte	0x4d
	.byte	0x25
	.long	0x3d0
	.uleb128 0x5
	.string	"var"
	.byte	0x4d
	.byte	0x2c
	.long	0x72
	.uleb128 0xa
	.string	"res"
	.byte	0x4f
	.byte	0x10
	.long	0x2e
	.byte	0
	.uleb128 0x12
	.long	.LASF59
	.byte	0x45
	.quad	.LFB58
	.quad	.LFE58-.LFB58
	.uleb128 0x1
	.byte	0x9c
	.long	0xc1f
	.uleb128 0x29
	.long	.LASF60
	.byte	0x45
	.byte	0x1f
	.long	0x3d5
	.uleb128 0x1
	.byte	0x55
	.uleb128 0x29
	.long	.LASF61
	.byte	0x45
	.byte	0x38
	.long	0x3d0
	.uleb128 0x1
	.byte	0x54
	.uleb128 0x25
	.string	"n"
	.byte	0x45
	.byte	0x4a
	.long	0xc2
	.long	.LLST26
	.long	.LVUS26
	.uleb128 0x3d
	.long	.LLRL27
	.long	0xc08
	.uleb128 0x1f
	.string	"i"
	.byte	0x47
	.byte	0xf
	.long	0xc2
	.long	.LLST28
	.long	.LVUS28
	.byte	0
	.uleb128 0x14
	.quad	.LVL42
	.uleb128 0x1
	.byte	0x30
	.uleb128 0x14
	.quad	.LVL44
	.uleb128 0x1
	.byte	0x30
	.byte	0
	.uleb128 0x12
	.long	.LASF62
	.byte	0x37
	.quad	.LFB57
	.quad	.LFE57-.LFB57
	.uleb128 0x1
	.byte	0x9c
	.long	0xc8e
	.uleb128 0x1d
	.long	.LASF61
	.byte	0x37
	.byte	0x22
	.long	0x3d0
	.long	.LLST23
	.long	.LVUS23
	.uleb128 0x25
	.string	"d"
	.byte	0x37
	.byte	0x3a
	.long	0x12c
	.long	.LLST24
	.long	.LVUS24
	.uleb128 0x1a
	.string	"n"
	.byte	0x37
	.byte	0x43
	.long	0xc2
	.uleb128 0x1
	.byte	0x51
	.uleb128 0x3e
	.quad	.LBB79
	.quad	.LBE79-.LBB79
	.uleb128 0x1f
	.string	"i"
	.byte	0x39
	.byte	0xf
	.long	0xc2
	.long	.LLST25
	.long	.LVUS25
	.byte	0
	.byte	0
	.uleb128 0x3f
	.long	.LASF63
	.byte	0x1
	.byte	0x25
	.byte	0x6
	.byte	0x1
	.long	0xcb9
	.uleb128 0x5
	.string	"x"
	.byte	0x25
	.byte	0x22
	.long	0x3d5
	.uleb128 0x5
	.string	"y"
	.byte	0x25
	.byte	0x2f
	.long	0x3d5
	.uleb128 0xa
	.string	"tmp"
	.byte	0x26
	.byte	0xe
	.long	0x3d0
	.byte	0
	.uleb128 0x12
	.long	.LASF64
	.byte	0x1d
	.quad	.LFB55
	.quad	.LFE55-.LFB55
	.uleb128 0x1
	.byte	0x9c
	.long	0xdc5
	.uleb128 0x1a
	.string	"a"
	.byte	0x1d
	.byte	0x1e
	.long	0xdc5
	.uleb128 0x1
	.byte	0x55
	.uleb128 0x1a
	.string	"b"
	.byte	0x1d
	.byte	0x2c
	.long	0xdc5
	.uleb128 0x1
	.byte	0x54
	.uleb128 0xa
	.string	"tmp"
	.byte	0x1e
	.byte	0xf
	.long	0x13d
	.uleb128 0x9
	.long	0xf2b
	.quad	.LBI73
	.byte	.LVU68
	.quad	.LBB73
	.quad	.LBE73-.LBB73
	.byte	0x1f
	.long	0xd39
	.uleb128 0x2
	.long	0xf52
	.long	.LLST15
	.long	.LVUS15
	.uleb128 0x2
	.long	0xf46
	.long	.LLST16
	.long	.LVUS16
	.uleb128 0x18
	.long	0xf3a
	.byte	0
	.uleb128 0x9
	.long	0xf2b
	.quad	.LBI75
	.byte	.LVU73
	.quad	.LBB75
	.quad	.LBE75-.LBB75
	.byte	0x20
	.long	0xd84
	.uleb128 0x2
	.long	0xf52
	.long	.LLST17
	.long	.LVUS17
	.uleb128 0x2
	.long	0xf46
	.long	.LLST18
	.long	.LVUS18
	.uleb128 0x2
	.long	0xf3a
	.long	.LLST19
	.long	.LVUS19
	.byte	0
	.uleb128 0x16
	.long	0xf2b
	.quad	.LBI77
	.byte	.LVU78
	.quad	.LBB77
	.quad	.LBE77-.LBB77
	.byte	0x21
	.byte	0x5
	.uleb128 0x2
	.long	0xf52
	.long	.LLST20
	.long	.LVUS20
	.uleb128 0x18
	.long	0xf46
	.uleb128 0x2
	.long	0xf3a
	.long	.LLST21
	.long	.LVUS21
	.byte	0
	.byte	0
	.uleb128 0xb
	.long	0x13d
	.uleb128 0x15
	.long	.LASF65
	.byte	0x13
	.byte	0xc
	.long	0x2e
	.byte	0x1
	.long	0xe0f
	.uleb128 0x5
	.string	"p1"
	.byte	0x13
	.byte	0x2b
	.long	0x12c
	.uleb128 0x5
	.string	"p2"
	.byte	0x13
	.byte	0x3b
	.long	0x12c
	.uleb128 0xa
	.string	"d"
	.byte	0x14
	.byte	0x10
	.long	0x2e
	.uleb128 0x21
	.uleb128 0xa
	.string	"i"
	.byte	0x15
	.byte	0x16
	.long	0x48
	.uleb128 0x21
	.uleb128 0xa
	.string	"dd"
	.byte	0x16
	.byte	0x14
	.long	0x2e
	.byte	0
	.byte	0
	.byte	0
	.uleb128 0x12
	.long	.LASF66
	.byte	0xb
	.quad	.LFB53
	.quad	.LFE53-.LFB53
	.uleb128 0x1
	.byte	0x9c
	.long	0xf2b
	.uleb128 0x1a
	.string	"a"
	.byte	0xb
	.byte	0xe
	.long	0x12c
	.uleb128 0x1
	.byte	0x55
	.uleb128 0x1a
	.string	"b"
	.byte	0xb
	.byte	0x14
	.long	0x12c
	.uleb128 0x1
	.byte	0x54
	.uleb128 0xa
	.string	"tmp"
	.byte	0xc
	.byte	0x7
	.long	0x2e
	.uleb128 0x9
	.long	0xf2b
	.quad	.LBI60
	.byte	.LVU4
	.quad	.LBB60
	.quad	.LBE60-.LBB60
	.byte	0xd
	.long	0xe97
	.uleb128 0x2
	.long	0xf52
	.long	.LLST0
	.long	.LVUS0
	.uleb128 0x2
	.long	0xf46
	.long	.LLST1
	.long	.LVUS1
	.uleb128 0x2
	.long	0xf3a
	.long	.LLST2
	.long	.LVUS2
	.byte	0
	.uleb128 0x9
	.long	0xf2b
	.quad	.LBI62
	.byte	.LVU9
	.quad	.LBB62
	.quad	.LBE62-.LBB62
	.byte	0xe
	.long	0xee2
	.uleb128 0x2
	.long	0xf52
	.long	.LLST3
	.long	.LVUS3
	.uleb128 0x2
	.long	0xf46
	.long	.LLST4
	.long	.LVUS4
	.uleb128 0x2
	.long	0xf3a
	.long	.LLST5
	.long	.LVUS5
	.byte	0
	.uleb128 0x16
	.long	0xf2b
	.quad	.LBI64
	.byte	.LVU14
	.quad	.LBB64
	.quad	.LBE64-.LBB64
	.byte	0xf
	.byte	0x5
	.uleb128 0x2
	.long	0xf52
	.long	.LLST6
	.long	.LVUS6
	.uleb128 0x2
	.long	0xf46
	.long	.LLST7
	.long	.LVUS7
	.uleb128 0x2
	.long	0xf3a
	.long	.LLST8
	.long	.LVUS8
	.byte	0
	.byte	0
	.uleb128 0x2a
	.long	.LASF70
	.byte	0x2
	.byte	0x1a
	.long	0x4f
	.long	0xf5f
	.uleb128 0x8
	.long	.LASF67
	.byte	0x2
	.byte	0x1a
	.byte	0x1
	.long	0x51
	.uleb128 0x8
	.long	.LASF68
	.byte	0x2
	.byte	0x1a
	.byte	0x1
	.long	0xb5
	.uleb128 0x8
	.long	.LASF69
	.byte	0x2
	.byte	0x1a
	.byte	0x1
	.long	0x35
	.byte	0
	.uleb128 0x2a
	.long	.LASF71
	.byte	0x3
	.byte	0x54
	.long	0x72
	.long	0xf7c
	.uleb128 0x8
	.long	.LASF72
	.byte	0x3
	.byte	0x54
	.byte	0x20
	.long	0x9d
	.uleb128 0x28
	.byte	0
	.uleb128 0x10
	.long	0xdca
	.quad	.LFB54
	.quad	.LFE54-.LFB54
	.uleb128 0x1
	.byte	0x9c
	.long	0x1001
	.uleb128 0x2
	.long	0xdda
	.long	.LLST9
	.long	.LVUS9
	.uleb128 0xf
	.long	0xde4
	.uleb128 0x1
	.byte	0x54
	.uleb128 0x3
	.long	0xdee
	.long	.LLST10
	.long	.LVUS10
	.uleb128 0x27
	.long	0xdf7
	.long	.LLRL11
	.long	0xfea
	.uleb128 0x3
	.long	0xdf8
	.long	.LLST12
	.long	.LVUS12
	.uleb128 0x26
	.long	0xe01
	.long	.LLRL13
	.uleb128 0x3
	.long	0xe02
	.long	.LLST14
	.long	.LVUS14
	.byte	0
	.byte	0
	.uleb128 0x14
	.quad	.LVL12
	.uleb128 0x1
	.byte	0x30
	.uleb128 0x14
	.quad	.LVL22
	.uleb128 0x1
	.byte	0x30
	.byte	0
	.uleb128 0x10
	.long	0xc8e
	.quad	.LFB56
	.quad	.LFE56-.LFB56
	.uleb128 0x1
	.byte	0x9c
	.long	0x1038
	.uleb128 0xf
	.long	0xc9b
	.uleb128 0x1
	.byte	0x55
	.uleb128 0xf
	.long	0xca4
	.uleb128 0x1
	.byte	0x54
	.uleb128 0x3
	.long	0xcad
	.long	.LLST22
	.long	.LVUS22
	.byte	0
	.uleb128 0x10
	.long	0xb6d
	.quad	.LFB59
	.quad	.LFE59-.LFB59
	.uleb128 0x1
	.byte	0x9c
	.long	0x1076
	.uleb128 0xf
	.long	0xb7d
	.uleb128 0x1
	.byte	0x55
	.uleb128 0xf
	.long	0xb86
	.uleb128 0x1
	.byte	0x54
	.uleb128 0xf
	.long	0xb8f
	.uleb128 0x1
	.byte	0x51
	.uleb128 0x3
	.long	0xb9a
	.long	.LLST29
	.long	.LVUS29
	.byte	0
	.uleb128 0x10
	.long	0x782
	.quad	.LFB61
	.quad	.LFE61-.LFB61
	.uleb128 0x1
	.byte	0x9c
	.long	0x11d6
	.uleb128 0xf
	.long	0x792
	.uleb128 0x1
	.byte	0x55
	.uleb128 0x2
	.long	0x79d
	.long	.LLST45
	.long	.LVUS45
	.uleb128 0x2
	.long	0x7a8
	.long	.LLST46
	.long	.LVUS46
	.uleb128 0x2
	.long	0x7b4
	.long	.LLST47
	.long	.LVUS47
	.uleb128 0x3
	.long	0x7c0
	.long	.LLST48
	.long	.LVUS48
	.uleb128 0x3
	.long	0x7cb
	.long	.LLST49
	.long	.LVUS49
	.uleb128 0x40
	.long	0x7d4
	.quad	.LBB127
	.quad	.LBE127-.LBB127
	.long	0x1193
	.uleb128 0x3
	.long	0x7d5
	.long	.LLST50
	.long	.LVUS50
	.uleb128 0x13
	.long	0xb6d
	.quad	.LBI128
	.byte	.LVU266
	.long	.LLRL51
	.byte	0x6a
	.byte	0xe
	.long	0x114a
	.uleb128 0x2
	.long	0xb8f
	.long	.LLST52
	.long	.LVUS52
	.uleb128 0x18
	.long	0xb86
	.uleb128 0x2
	.long	0xb7d
	.long	.LLST53
	.long	.LVUS53
	.uleb128 0xe
	.long	.LLRL51
	.uleb128 0x3
	.long	0xb9a
	.long	.LLST54
	.long	.LVUS54
	.byte	0
	.byte	0
	.uleb128 0x16
	.long	0xc8e
	.quad	.LBI134
	.byte	.LVU277
	.quad	.LBB134
	.quad	.LBE134-.LBB134
	.byte	0x6c
	.byte	0xd
	.uleb128 0x2
	.long	0xca4
	.long	.LLST55
	.long	.LVUS55
	.uleb128 0x2
	.long	0xc9b
	.long	.LLST56
	.long	.LVUS56
	.uleb128 0x3
	.long	0xcad
	.long	.LLST57
	.long	.LVUS57
	.byte	0
	.byte	0
	.uleb128 0x19
	.long	0xc8e
	.quad	.LBI136
	.byte	.LVU293
	.long	.LLRL58
	.byte	0x6f
	.byte	0x5
	.uleb128 0x2
	.long	0xca4
	.long	.LLST59
	.long	.LVUS59
	.uleb128 0x2
	.long	0xc9b
	.long	.LLST60
	.long	.LVUS60
	.uleb128 0xe
	.long	.LLRL58
	.uleb128 0x3
	.long	0xcad
	.long	.LLST61
	.long	.LVUS61
	.byte	0
	.byte	0
	.byte	0
	.uleb128 0x10
	.long	0x72e
	.quad	.LFB62
	.quad	.LFE62-.LFB62
	.uleb128 0x1
	.byte	0x9c
	.long	0x149c
	.uleb128 0x2
	.long	0x73e
	.long	.LLST62
	.long	.LVUS62
	.uleb128 0x2
	.long	0x747
	.long	.LLST63
	.long	.LVUS63
	.uleb128 0x2
	.long	0x753
	.long	.LLST64
	.long	.LVUS64
	.uleb128 0x2
	.long	0x75f
	.long	.LLST65
	.long	.LVUS65
	.uleb128 0x3
	.long	0x76b
	.long	.LLST66
	.long	.LVUS66
	.uleb128 0x27
	.long	0x774
	.long	.LLRL67
	.long	0x13a6
	.uleb128 0x3
	.long	0x775
	.long	.LLST68
	.long	.LVUS68
	.uleb128 0x19
	.long	0x782
	.quad	.LBI158
	.byte	.LVU327
	.long	.LLRL69
	.byte	0x8f
	.byte	0x1a
	.uleb128 0x2
	.long	0x7b4
	.long	.LLST70
	.long	.LVUS70
	.uleb128 0x2
	.long	0x7a8
	.long	.LLST71
	.long	.LVUS71
	.uleb128 0x2
	.long	0x79d
	.long	.LLST72
	.long	.LVUS72
	.uleb128 0x2
	.long	0x792
	.long	.LLST73
	.long	.LVUS73
	.uleb128 0xe
	.long	.LLRL69
	.uleb128 0x3
	.long	0x7c0
	.long	.LLST74
	.long	.LVUS74
	.uleb128 0x3
	.long	0x7cb
	.long	.LLST75
	.long	.LVUS75
	.uleb128 0x27
	.long	0x7d4
	.long	.LLRL76
	.long	0x1361
	.uleb128 0x3
	.long	0x7d5
	.long	.LLST77
	.long	.LVUS77
	.uleb128 0x13
	.long	0xb6d
	.quad	.LBI161
	.byte	.LVU342
	.long	.LLRL78
	.byte	0x6a
	.byte	0xe
	.long	0x1318
	.uleb128 0x2
	.long	0xb8f
	.long	.LLST79
	.long	.LVUS79
	.uleb128 0x18
	.long	0xb86
	.uleb128 0x2
	.long	0xb7d
	.long	.LLST80
	.long	.LVUS80
	.uleb128 0xe
	.long	.LLRL78
	.uleb128 0x3
	.long	0xb9a
	.long	.LLST81
	.long	.LVUS81
	.byte	0
	.byte	0
	.uleb128 0x16
	.long	0xc8e
	.quad	.LBI169
	.byte	.LVU353
	.quad	.LBB169
	.quad	.LBE169-.LBB169
	.byte	0x6c
	.byte	0xd
	.uleb128 0x2
	.long	0xca4
	.long	.LLST82
	.long	.LVUS82
	.uleb128 0x2
	.long	0xc9b
	.long	.LLST83
	.long	.LVUS83
	.uleb128 0x3
	.long	0xcad
	.long	.LLST84
	.long	.LVUS84
	.byte	0
	.byte	0
	.uleb128 0x19
	.long	0xc8e
	.quad	.LBI173
	.byte	.LVU370
	.long	.LLRL85
	.byte	0x6f
	.byte	0x5
	.uleb128 0x2
	.long	0xca4
	.long	.LLST86
	.long	.LVUS86
	.uleb128 0x2
	.long	0xc9b
	.long	.LLST87
	.long	.LVUS87
	.uleb128 0xe
	.long	.LLRL85
	.uleb128 0x3
	.long	0xcad
	.long	.LLST88
	.long	.LVUS88
	.byte	0
	.byte	0
	.byte	0
	.byte	0
	.byte	0
	.uleb128 0x16
	.long	0x72e
	.quad	.LBI187
	.byte	.LVU405
	.quad	.LBB187
	.quad	.LBE187-.LBB187
	.byte	0x74
	.byte	0x5
	.uleb128 0x2
	.long	0x75f
	.long	.LLST89
	.long	.LVUS89
	.uleb128 0x2
	.long	0x753
	.long	.LLST90
	.long	.LVUS90
	.uleb128 0x2
	.long	0x747
	.long	.LLST91
	.long	.LVUS91
	.uleb128 0x2
	.long	0x73e
	.long	.LLST92
	.long	.LVUS92
	.uleb128 0x41
	.long	0x76b
	.uleb128 0x13
	.long	0xb6d
	.quad	.LBI189
	.byte	.LVU416
	.long	.LLRL93
	.byte	0x83
	.byte	0xc
	.long	0x1452
	.uleb128 0x2
	.long	0xb8f
	.long	.LLST94
	.long	.LVUS94
	.uleb128 0x2
	.long	0xb86
	.long	.LLST95
	.long	.LVUS95
	.uleb128 0x2
	.long	0xb7d
	.long	.LLST96
	.long	.LVUS96
	.uleb128 0xe
	.long	.LLRL93
	.uleb128 0x3
	.long	0xb9a
	.long	.LLST97
	.long	.LVUS97
	.byte	0
	.byte	0
	.uleb128 0x16
	.long	0xc8e
	.quad	.LBI193
	.byte	.LVU424
	.quad	.LBB193
	.quad	.LBE193-.LBB193
	.byte	0x83
	.byte	0x35
	.uleb128 0x2
	.long	0xca4
	.long	.LLST98
	.long	.LVUS98
	.uleb128 0x2
	.long	0xc9b
	.long	.LLST99
	.long	.LVUS99
	.uleb128 0x3
	.long	0xcad
	.long	.LLST100
	.long	.LVUS100
	.byte	0
	.byte	0
	.byte	0
	.uleb128 0x10
	.long	0x6c6
	.quad	.LFB63
	.quad	.LFE63-.LFB63
	.uleb128 0x1
	.byte	0x9c
	.long	0x15b8
	.uleb128 0x2
	.long	0x6d6
	.long	.LLST101
	.long	.LVUS101
	.uleb128 0x2
	.long	0x6df
	.long	.LLST102
	.long	.LVUS102
	.uleb128 0x2
	.long	0x6eb
	.long	.LLST103
	.long	.LVUS103
	.uleb128 0x2
	.long	0x6f6
	.long	.LLST104
	.long	.LVUS104
	.uleb128 0x2
	.long	0x702
	.long	.LLST105
	.long	.LVUS105
	.uleb128 0x3
	.long	0x70e
	.long	.LLST106
	.long	.LVUS106
	.uleb128 0x3
	.long	0x717
	.long	.LLST107
	.long	.LVUS107
	.uleb128 0x3
	.long	0x722
	.long	.LLST108
	.long	.LVUS108
	.uleb128 0x4
	.quad	.LVL124
	.long	0x72e
	.long	0x154e
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x55
	.uleb128 0x4
	.byte	0x91
	.sleb128 -80
	.byte	0x6
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x54
	.uleb128 0x5
	.byte	0x91
	.sleb128 -72
	.byte	0x94
	.byte	0x4
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x51
	.uleb128 0x2
	.byte	0x7f
	.sleb128 0
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x52
	.uleb128 0x2
	.byte	0x7e
	.sleb128 0
	.byte	0
	.uleb128 0x4
	.quad	.LVL128
	.long	0x6c6
	.long	0x1586
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x55
	.uleb128 0x4
	.byte	0x91
	.sleb128 -80
	.byte	0x6
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x54
	.uleb128 0x5
	.byte	0x91
	.sleb128 -72
	.byte	0x94
	.byte	0x4
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x51
	.uleb128 0x2
	.byte	0x7c
	.sleb128 -1
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x52
	.uleb128 0x2
	.byte	0x73
	.sleb128 0
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x58
	.uleb128 0x5
	.byte	0x91
	.sleb128 -68
	.byte	0x94
	.byte	0x4
	.byte	0
	.uleb128 0x6
	.quad	.LVL129
	.long	0x6c6
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x55
	.uleb128 0x4
	.byte	0x91
	.sleb128 -80
	.byte	0x6
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x54
	.uleb128 0x2
	.byte	0x7c
	.sleb128 1
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x51
	.uleb128 0x2
	.byte	0x7f
	.sleb128 0
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x52
	.uleb128 0x2
	.byte	0x73
	.sleb128 0
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x58
	.uleb128 0x5
	.byte	0x91
	.sleb128 -68
	.byte	0x94
	.byte	0x4
	.byte	0
	.byte	0
	.uleb128 0x10
	.long	0x6c6
	.quad	.LFB71
	.quad	.LFE71-.LFB71
	.uleb128 0x1
	.byte	0x9c
	.long	0x16c3
	.uleb128 0x2
	.long	0x6d6
	.long	.LLST109
	.long	.LVUS109
	.uleb128 0x2
	.long	0x6df
	.long	.LLST110
	.long	.LVUS110
	.uleb128 0x2
	.long	0x6eb
	.long	.LLST111
	.long	.LVUS111
	.uleb128 0x2
	.long	0x6f6
	.long	.LLST112
	.long	.LVUS112
	.uleb128 0x3
	.long	0x70e
	.long	.LLST113
	.long	.LVUS113
	.uleb128 0x3
	.long	0x717
	.long	.LLST114
	.long	.LVUS114
	.uleb128 0x3
	.long	0x722
	.long	.LLST115
	.long	.LVUS115
	.uleb128 0x2b
	.long	0x702
	.byte	0x2
	.uleb128 0x4
	.quad	.LVL140
	.long	0x72e
	.long	0x165e
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x55
	.uleb128 0x2
	.byte	0x7c
	.sleb128 0
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x54
	.uleb128 0x2
	.byte	0x7d
	.sleb128 0
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x51
	.uleb128 0x2
	.byte	0x76
	.sleb128 0
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x52
	.uleb128 0x2
	.byte	0x7f
	.sleb128 0
	.byte	0
	.uleb128 0x4
	.quad	.LVL144
	.long	0x6c6
	.long	0x1692
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x55
	.uleb128 0x2
	.byte	0x7c
	.sleb128 0
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x54
	.uleb128 0x2
	.byte	0x7d
	.sleb128 0
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x51
	.uleb128 0x7
	.byte	0x91
	.sleb128 -68
	.byte	0x94
	.byte	0x4
	.byte	0x31
	.byte	0x1c
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x52
	.uleb128 0x2
	.byte	0x73
	.sleb128 0
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x58
	.uleb128 0x1
	.byte	0x33
	.byte	0
	.uleb128 0x6
	.quad	.LVL145
	.long	0x6c6
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x55
	.uleb128 0x2
	.byte	0x7c
	.sleb128 0
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x54
	.uleb128 0x7
	.byte	0x91
	.sleb128 -68
	.byte	0x94
	.byte	0x4
	.byte	0x23
	.uleb128 0x1
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x51
	.uleb128 0x2
	.byte	0x76
	.sleb128 0
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x52
	.uleb128 0x2
	.byte	0x73
	.sleb128 0
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x58
	.uleb128 0x1
	.byte	0x33
	.byte	0
	.byte	0
	.uleb128 0x10
	.long	0x6c6
	.quad	.LFB72
	.quad	.LFE72-.LFB72
	.uleb128 0x1
	.byte	0x9c
	.long	0x17d2
	.uleb128 0x2
	.long	0x6d6
	.long	.LLST116
	.long	.LVUS116
	.uleb128 0x2
	.long	0x6df
	.long	.LLST117
	.long	.LVUS117
	.uleb128 0x2
	.long	0x6eb
	.long	.LLST118
	.long	.LVUS118
	.uleb128 0x2
	.long	0x6f6
	.long	.LLST119
	.long	.LVUS119
	.uleb128 0x3
	.long	0x70e
	.long	.LLST120
	.long	.LVUS120
	.uleb128 0x3
	.long	0x717
	.long	.LLST121
	.long	.LVUS121
	.uleb128 0x3
	.long	0x722
	.long	.LLST122
	.long	.LVUS122
	.uleb128 0x2b
	.long	0x702
	.byte	0x1
	.uleb128 0x4
	.quad	.LVL156
	.long	0x72e
	.long	0x1769
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x55
	.uleb128 0x2
	.byte	0x7c
	.sleb128 0
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x54
	.uleb128 0x2
	.byte	0x7d
	.sleb128 0
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x51
	.uleb128 0x2
	.byte	0x76
	.sleb128 0
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x52
	.uleb128 0x2
	.byte	0x7f
	.sleb128 0
	.byte	0
	.uleb128 0x4
	.quad	.LVL160
	.long	0x15b8
	.long	0x179f
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x55
	.uleb128 0x2
	.byte	0x7c
	.sleb128 0
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x54
	.uleb128 0x2
	.byte	0x7d
	.sleb128 0
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x51
	.uleb128 0x7
	.byte	0x91
	.sleb128 -68
	.byte	0x94
	.byte	0x4
	.byte	0x31
	.byte	0x1c
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x52
	.uleb128 0x2
	.byte	0x73
	.sleb128 0
	.uleb128 0x1c
	.long	0x702
	.uleb128 0x1
	.byte	0x32
	.byte	0
	.uleb128 0x6
	.quad	.LVL161
	.long	0x15b8
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x55
	.uleb128 0x2
	.byte	0x7c
	.sleb128 0
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x54
	.uleb128 0x7
	.byte	0x91
	.sleb128 -68
	.byte	0x94
	.byte	0x4
	.byte	0x23
	.uleb128 0x1
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x51
	.uleb128 0x2
	.byte	0x76
	.sleb128 0
	.uleb128 0x1
	.uleb128 0x1
	.byte	0x52
	.uleb128 0x2
	.byte	0x73
	.sleb128 0
	.uleb128 0x1c
	.long	0x702
	.uleb128 0x1
	.byte	0x32
	.byte	0
	.byte	0
	.uleb128 0x42
	.long	.LASF74
	.long	.LASF74
	.uleb128 0x43
	.uleb128 0x4
	.byte	0x9e
	.uleb128 0x2
	.byte	0xa
	.byte	0
	.uleb128 0x44
	.long	.LASF75
	.long	.LASF76
	.byte	0xa
	.byte	0
	.byte	0
	.section	.debug_abbrev,"",@progbits
.Ldebug_abbrev0:
	.uleb128 0x1
	.uleb128 0x49
	.byte	0
	.uleb128 0x2
	.uleb128 0x18
	.uleb128 0x7e
	.uleb128 0x18
	.byte	0
	.byte	0
	.uleb128 0x2
	.uleb128 0x5
	.byte	0
	.uleb128 0x31
	.uleb128 0x13
	.uleb128 0x2
	.uleb128 0x17
	.uleb128 0x2137
	.uleb128 0x17
	.byte	0
	.byte	0
	.uleb128 0x3
	.uleb128 0x34
	.byte	0
	.uleb128 0x31
	.uleb128 0x13
	.uleb128 0x2
	.uleb128 0x17
	.uleb128 0x2137
	.uleb128 0x17
	.byte	0
	.byte	0
	.uleb128 0x4
	.uleb128 0x48
	.byte	0x1
	.uleb128 0x7d
	.uleb128 0x1
	.uleb128 0x7f
	.uleb128 0x13
	.uleb128 0x1
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0x5
	.uleb128 0x5
	.byte	0
	.uleb128 0x3
	.uleb128 0x8
	.uleb128 0x3a
	.uleb128 0x21
	.sleb128 1
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x39
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0x6
	.uleb128 0x48
	.byte	0x1
	.uleb128 0x7d
	.uleb128 0x1
	.uleb128 0x7f
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0x7
	.uleb128 0x24
	.byte	0
	.uleb128 0xb
	.uleb128 0xb
	.uleb128 0x3e
	.uleb128 0xb
	.uleb128 0x3
	.uleb128 0xe
	.byte	0
	.byte	0
	.uleb128 0x8
	.uleb128 0x5
	.byte	0
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x39
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0x9
	.uleb128 0x1d
	.byte	0x1
	.uleb128 0x31
	.uleb128 0x13
	.uleb128 0x52
	.uleb128 0x1
	.uleb128 0x2138
	.uleb128 0xb
	.uleb128 0x11
	.uleb128 0x1
	.uleb128 0x12
	.uleb128 0x7
	.uleb128 0x58
	.uleb128 0x21
	.sleb128 1
	.uleb128 0x59
	.uleb128 0xb
	.uleb128 0x57
	.uleb128 0x21
	.sleb128 5
	.uleb128 0x1
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0xa
	.uleb128 0x34
	.byte	0
	.uleb128 0x3
	.uleb128 0x8
	.uleb128 0x3a
	.uleb128 0x21
	.sleb128 1
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x39
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0xb
	.uleb128 0xf
	.byte	0
	.uleb128 0xb
	.uleb128 0x21
	.sleb128 8
	.uleb128 0x49
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0xc
	.uleb128 0xd
	.byte	0
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x39
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x38
	.uleb128 0xb
	.byte	0
	.byte	0
	.uleb128 0xd
	.uleb128 0x5
	.byte	0
	.uleb128 0x49
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0xe
	.uleb128 0xb
	.byte	0x1
	.uleb128 0x55
	.uleb128 0x17
	.byte	0
	.byte	0
	.uleb128 0xf
	.uleb128 0x5
	.byte	0
	.uleb128 0x31
	.uleb128 0x13
	.uleb128 0x2
	.uleb128 0x18
	.byte	0
	.byte	0
	.uleb128 0x10
	.uleb128 0x2e
	.byte	0x1
	.uleb128 0x31
	.uleb128 0x13
	.uleb128 0x11
	.uleb128 0x1
	.uleb128 0x12
	.uleb128 0x7
	.uleb128 0x40
	.uleb128 0x18
	.uleb128 0x7a
	.uleb128 0x19
	.uleb128 0x1
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0x11
	.uleb128 0x16
	.byte	0
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x39
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0x12
	.uleb128 0x2e
	.byte	0x1
	.uleb128 0x3f
	.uleb128 0x19
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0x21
	.sleb128 1
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x39
	.uleb128 0x21
	.sleb128 6
	.uleb128 0x27
	.uleb128 0x19
	.uleb128 0x11
	.uleb128 0x1
	.uleb128 0x12
	.uleb128 0x7
	.uleb128 0x40
	.uleb128 0x18
	.uleb128 0x7a
	.uleb128 0x19
	.uleb128 0x1
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0x13
	.uleb128 0x1d
	.byte	0x1
	.uleb128 0x31
	.uleb128 0x13
	.uleb128 0x52
	.uleb128 0x1
	.uleb128 0x2138
	.uleb128 0xb
	.uleb128 0x55
	.uleb128 0x17
	.uleb128 0x58
	.uleb128 0x21
	.sleb128 1
	.uleb128 0x59
	.uleb128 0xb
	.uleb128 0x57
	.uleb128 0xb
	.uleb128 0x1
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0x14
	.uleb128 0x48
	.byte	0
	.uleb128 0x7d
	.uleb128 0x1
	.uleb128 0x83
	.uleb128 0x18
	.byte	0
	.byte	0
	.uleb128 0x15
	.uleb128 0x2e
	.byte	0x1
	.uleb128 0x3f
	.uleb128 0x19
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0x21
	.sleb128 1
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x39
	.uleb128 0xb
	.uleb128 0x27
	.uleb128 0x19
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x20
	.uleb128 0xb
	.uleb128 0x1
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0x16
	.uleb128 0x1d
	.byte	0x1
	.uleb128 0x31
	.uleb128 0x13
	.uleb128 0x52
	.uleb128 0x1
	.uleb128 0x2138
	.uleb128 0xb
	.uleb128 0x11
	.uleb128 0x1
	.uleb128 0x12
	.uleb128 0x7
	.uleb128 0x58
	.uleb128 0x21
	.sleb128 1
	.uleb128 0x59
	.uleb128 0xb
	.uleb128 0x57
	.uleb128 0xb
	.byte	0
	.byte	0
	.uleb128 0x17
	.uleb128 0x5
	.byte	0
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0x21
	.sleb128 1
	.uleb128 0x3b
	.uleb128 0x5
	.uleb128 0x39
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x2
	.uleb128 0x17
	.uleb128 0x2137
	.uleb128 0x17
	.byte	0
	.byte	0
	.uleb128 0x18
	.uleb128 0x5
	.byte	0
	.uleb128 0x31
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0x19
	.uleb128 0x1d
	.byte	0x1
	.uleb128 0x31
	.uleb128 0x13
	.uleb128 0x52
	.uleb128 0x1
	.uleb128 0x2138
	.uleb128 0xb
	.uleb128 0x55
	.uleb128 0x17
	.uleb128 0x58
	.uleb128 0x21
	.sleb128 1
	.uleb128 0x59
	.uleb128 0xb
	.uleb128 0x57
	.uleb128 0xb
	.byte	0
	.byte	0
	.uleb128 0x1a
	.uleb128 0x5
	.byte	0
	.uleb128 0x3
	.uleb128 0x8
	.uleb128 0x3a
	.uleb128 0x21
	.sleb128 1
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x39
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x2
	.uleb128 0x18
	.byte	0
	.byte	0
	.uleb128 0x1b
	.uleb128 0x2e
	.byte	0x1
	.uleb128 0x3f
	.uleb128 0x19
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0x21
	.sleb128 7
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x39
	.uleb128 0x21
	.sleb128 6
	.uleb128 0x27
	.uleb128 0x19
	.uleb128 0x3c
	.uleb128 0x19
	.uleb128 0x1
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0x1c
	.uleb128 0x49
	.byte	0
	.uleb128 0x80
	.uleb128 0x13
	.uleb128 0x7e
	.uleb128 0x18
	.byte	0
	.byte	0
	.uleb128 0x1d
	.uleb128 0x5
	.byte	0
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0x21
	.sleb128 1
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x39
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x2
	.uleb128 0x17
	.uleb128 0x2137
	.uleb128 0x17
	.byte	0
	.byte	0
	.uleb128 0x1e
	.uleb128 0x34
	.byte	0
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0x21
	.sleb128 1
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x39
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x2
	.uleb128 0x17
	.uleb128 0x2137
	.uleb128 0x17
	.byte	0
	.byte	0
	.uleb128 0x1f
	.uleb128 0x34
	.byte	0
	.uleb128 0x3
	.uleb128 0x8
	.uleb128 0x3a
	.uleb128 0x21
	.sleb128 1
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x39
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x2
	.uleb128 0x17
	.uleb128 0x2137
	.uleb128 0x17
	.byte	0
	.byte	0
	.uleb128 0x20
	.uleb128 0x34
	.byte	0
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0x21
	.sleb128 1
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x39
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0x21
	.uleb128 0xb
	.byte	0x1
	.byte	0
	.byte	0
	.uleb128 0x22
	.uleb128 0x37
	.byte	0
	.uleb128 0x49
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0x23
	.uleb128 0x13
	.byte	0x1
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0xb
	.uleb128 0xb
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x39
	.uleb128 0x21
	.sleb128 8
	.uleb128 0x1
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0x24
	.uleb128 0xd
	.byte	0
	.uleb128 0x3
	.uleb128 0x8
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x39
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x38
	.uleb128 0xb
	.byte	0
	.byte	0
	.uleb128 0x25
	.uleb128 0x5
	.byte	0
	.uleb128 0x3
	.uleb128 0x8
	.uleb128 0x3a
	.uleb128 0x21
	.sleb128 1
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x39
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x2
	.uleb128 0x17
	.uleb128 0x2137
	.uleb128 0x17
	.byte	0
	.byte	0
	.uleb128 0x26
	.uleb128 0xb
	.byte	0x1
	.uleb128 0x31
	.uleb128 0x13
	.uleb128 0x55
	.uleb128 0x17
	.byte	0
	.byte	0
	.uleb128 0x27
	.uleb128 0xb
	.byte	0x1
	.uleb128 0x31
	.uleb128 0x13
	.uleb128 0x55
	.uleb128 0x17
	.uleb128 0x1
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0x28
	.uleb128 0x18
	.byte	0
	.byte	0
	.byte	0
	.uleb128 0x29
	.uleb128 0x5
	.byte	0
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0x21
	.sleb128 1
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x39
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x2
	.uleb128 0x18
	.byte	0
	.byte	0
	.uleb128 0x2a
	.uleb128 0x2e
	.byte	0x1
	.uleb128 0x3f
	.uleb128 0x19
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x39
	.uleb128 0x21
	.sleb128 1
	.uleb128 0x27
	.uleb128 0x19
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x20
	.uleb128 0x21
	.sleb128 3
	.uleb128 0x34
	.uleb128 0x19
	.uleb128 0x1
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0x2b
	.uleb128 0x5
	.byte	0
	.uleb128 0x31
	.uleb128 0x13
	.uleb128 0x1c
	.uleb128 0xb
	.byte	0
	.byte	0
	.uleb128 0x2c
	.uleb128 0x11
	.byte	0x1
	.uleb128 0x25
	.uleb128 0xe
	.uleb128 0x13
	.uleb128 0xb
	.uleb128 0x3
	.uleb128 0x1f
	.uleb128 0x1b
	.uleb128 0x1f
	.uleb128 0x11
	.uleb128 0x1
	.uleb128 0x12
	.uleb128 0x7
	.uleb128 0x10
	.uleb128 0x17
	.byte	0
	.byte	0
	.uleb128 0x2d
	.uleb128 0xf
	.byte	0
	.uleb128 0xb
	.uleb128 0xb
	.byte	0
	.byte	0
	.uleb128 0x2e
	.uleb128 0x24
	.byte	0
	.uleb128 0xb
	.uleb128 0xb
	.uleb128 0x3e
	.uleb128 0xb
	.uleb128 0x3
	.uleb128 0x8
	.byte	0
	.byte	0
	.uleb128 0x2f
	.uleb128 0x26
	.byte	0
	.uleb128 0x49
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0x30
	.uleb128 0x26
	.byte	0
	.byte	0
	.byte	0
	.uleb128 0x31
	.uleb128 0x34
	.byte	0
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x39
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x3f
	.uleb128 0x19
	.uleb128 0x3c
	.uleb128 0x19
	.byte	0
	.byte	0
	.uleb128 0x32
	.uleb128 0x34
	.byte	0
	.uleb128 0x3
	.uleb128 0x8
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x39
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x3f
	.uleb128 0x19
	.uleb128 0x2
	.uleb128 0x18
	.byte	0
	.byte	0
	.uleb128 0x33
	.uleb128 0x2e
	.byte	0x1
	.uleb128 0x3f
	.uleb128 0x19
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x39
	.uleb128 0xb
	.uleb128 0x27
	.uleb128 0x19
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x3c
	.uleb128 0x19
	.uleb128 0x1
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0x34
	.uleb128 0x2e
	.byte	0x1
	.uleb128 0x3f
	.uleb128 0x19
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0x5
	.uleb128 0x39
	.uleb128 0xb
	.uleb128 0x27
	.uleb128 0x19
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x11
	.uleb128 0x1
	.uleb128 0x12
	.uleb128 0x7
	.uleb128 0x40
	.uleb128 0x18
	.uleb128 0x7a
	.uleb128 0x19
	.uleb128 0x1
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0x35
	.uleb128 0x5
	.byte	0
	.uleb128 0x3
	.uleb128 0x8
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0x5
	.uleb128 0x39
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x2
	.uleb128 0x17
	.uleb128 0x2137
	.uleb128 0x17
	.byte	0
	.byte	0
	.uleb128 0x36
	.uleb128 0x34
	.byte	0
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0x5
	.uleb128 0x39
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x2
	.uleb128 0x17
	.uleb128 0x2137
	.uleb128 0x17
	.byte	0
	.byte	0
	.uleb128 0x37
	.uleb128 0x1d
	.byte	0x1
	.uleb128 0x31
	.uleb128 0x13
	.uleb128 0x52
	.uleb128 0x1
	.uleb128 0x2138
	.uleb128 0xb
	.uleb128 0x55
	.uleb128 0x17
	.uleb128 0x58
	.uleb128 0xb
	.uleb128 0x59
	.uleb128 0x5
	.uleb128 0x57
	.uleb128 0xb
	.byte	0
	.byte	0
	.uleb128 0x38
	.uleb128 0x2e
	.byte	0x1
	.uleb128 0x3f
	.uleb128 0x19
	.uleb128 0x3
	.uleb128 0x8
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0x5
	.uleb128 0x39
	.uleb128 0xb
	.uleb128 0x27
	.uleb128 0x19
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x11
	.uleb128 0x1
	.uleb128 0x12
	.uleb128 0x7
	.uleb128 0x40
	.uleb128 0x18
	.uleb128 0x7a
	.uleb128 0x19
	.uleb128 0x1
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0x39
	.uleb128 0x34
	.byte	0
	.uleb128 0x3
	.uleb128 0x8
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0x5
	.uleb128 0x39
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x2
	.uleb128 0x18
	.byte	0
	.byte	0
	.uleb128 0x3a
	.uleb128 0x48
	.byte	0
	.uleb128 0x7d
	.uleb128 0x1
	.uleb128 0x7f
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0x3b
	.uleb128 0xb
	.byte	0x1
	.uleb128 0x11
	.uleb128 0x1
	.uleb128 0x12
	.uleb128 0x7
	.uleb128 0x1
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0x3c
	.uleb128 0x48
	.byte	0x1
	.uleb128 0x7d
	.uleb128 0x1
	.uleb128 0x82
	.uleb128 0x19
	.uleb128 0x7f
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0x3d
	.uleb128 0xb
	.byte	0x1
	.uleb128 0x55
	.uleb128 0x17
	.uleb128 0x1
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0x3e
	.uleb128 0xb
	.byte	0x1
	.uleb128 0x11
	.uleb128 0x1
	.uleb128 0x12
	.uleb128 0x7
	.byte	0
	.byte	0
	.uleb128 0x3f
	.uleb128 0x2e
	.byte	0x1
	.uleb128 0x3f
	.uleb128 0x19
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x39
	.uleb128 0xb
	.uleb128 0x27
	.uleb128 0x19
	.uleb128 0x20
	.uleb128 0xb
	.uleb128 0x1
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0x40
	.uleb128 0xb
	.byte	0x1
	.uleb128 0x31
	.uleb128 0x13
	.uleb128 0x11
	.uleb128 0x1
	.uleb128 0x12
	.uleb128 0x7
	.uleb128 0x1
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0x41
	.uleb128 0x34
	.byte	0
	.uleb128 0x31
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0x42
	.uleb128 0x2e
	.byte	0
	.uleb128 0x3f
	.uleb128 0x19
	.uleb128 0x3c
	.uleb128 0x19
	.uleb128 0x6e
	.uleb128 0xe
	.uleb128 0x3
	.uleb128 0xe
	.byte	0
	.byte	0
	.uleb128 0x43
	.uleb128 0x36
	.byte	0
	.uleb128 0x2
	.uleb128 0x18
	.byte	0
	.byte	0
	.uleb128 0x44
	.uleb128 0x2e
	.byte	0
	.uleb128 0x3f
	.uleb128 0x19
	.uleb128 0x3c
	.uleb128 0x19
	.uleb128 0x6e
	.uleb128 0xe
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.byte	0
	.byte	0
	.byte	0
	.section	.debug_loclists,"",@progbits
	.long	.Ldebug_loc3-.Ldebug_loc2
.Ldebug_loc2:
	.value	0x5
	.byte	0x8
	.byte	0
	.long	0
.Ldebug_loc0:
.LVUS145:
	.uleb128 0
	.uleb128 .LVU758
	.uleb128 .LVU758
	.uleb128 .LVU782
	.uleb128 .LVU782
	.uleb128 .LVU787
	.uleb128 .LVU787
	.uleb128 0
.LLST145:
	.byte	0x4
	.uleb128 .LVL216-.Ltext0
	.uleb128 .LVL221-1-.Ltext0
	.uleb128 0x1
	.byte	0x55
	.byte	0x4
	.uleb128 .LVL221-1-.Ltext0
	.uleb128 .LVL229-.Ltext0
	.uleb128 0x1
	.byte	0x56
	.byte	0x4
	.uleb128 .LVL229-.Ltext0
	.uleb128 .LVL232-.Ltext0
	.uleb128 0x4
	.byte	0xa3
	.uleb128 0x1
	.byte	0x55
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL232-.Ltext0
	.uleb128 .LFE68-.Ltext0
	.uleb128 0x1
	.byte	0x55
	.byte	0
.LVUS146:
	.uleb128 0
	.uleb128 .LVU757
	.uleb128 .LVU757
	.uleb128 .LVU787
	.uleb128 .LVU787
	.uleb128 0
.LLST146:
	.byte	0x4
	.uleb128 .LVL216-.Ltext0
	.uleb128 .LVL220-.Ltext0
	.uleb128 0x1
	.byte	0x54
	.byte	0x4
	.uleb128 .LVL220-.Ltext0
	.uleb128 .LVL232-.Ltext0
	.uleb128 0x4
	.byte	0xa3
	.uleb128 0x1
	.byte	0x54
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL232-.Ltext0
	.uleb128 .LFE68-.Ltext0
	.uleb128 0x1
	.byte	0x54
	.byte	0
.LVUS147:
	.uleb128 0
	.uleb128 .LVU756
	.uleb128 .LVU756
	.uleb128 .LVU787
	.uleb128 .LVU787
	.uleb128 0
.LLST147:
	.byte	0x4
	.uleb128 .LVL216-.Ltext0
	.uleb128 .LVL219-.Ltext0
	.uleb128 0x1
	.byte	0x51
	.byte	0x4
	.uleb128 .LVL219-.Ltext0
	.uleb128 .LVL232-.Ltext0
	.uleb128 0x4
	.byte	0xa3
	.uleb128 0x1
	.byte	0x51
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL232-.Ltext0
	.uleb128 .LFE68-.Ltext0
	.uleb128 0x1
	.byte	0x51
	.byte	0
.LVUS148:
	.uleb128 .LVU778
	.uleb128 .LVU781
	.uleb128 .LVU781
	.uleb128 .LVU783
.LLST148:
	.byte	0x4
	.uleb128 .LVL227-.Ltext0
	.uleb128 .LVL228-.Ltext0
	.uleb128 0x1
	.byte	0x53
	.byte	0x4
	.uleb128 .LVL228-.Ltext0
	.uleb128 .LVL230-.Ltext0
	.uleb128 0x1
	.byte	0x50
	.byte	0
.LVUS150:
	.uleb128 .LVU746
	.uleb128 .LVU778
	.uleb128 .LVU783
	.uleb128 0
.LLST150:
	.byte	0x4
	.uleb128 .LVL217-.Ltext0
	.uleb128 .LVL227-.Ltext0
	.uleb128 0x2
	.byte	0x30
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL230-.Ltext0
	.uleb128 .LFE68-.Ltext0
	.uleb128 0x2
	.byte	0x30
	.byte	0x9f
	.byte	0
.LVUS153:
	.uleb128 .LVU745
	.uleb128 .LVU751
	.uleb128 .LVU751
	.uleb128 .LVU778
	.uleb128 .LVU783
	.uleb128 .LVU786
	.uleb128 .LVU786
	.uleb128 .LVU787
	.uleb128 .LVU787
	.uleb128 .LVU802
	.uleb128 .LVU802
	.uleb128 0
.LLST153:
	.byte	0x4
	.uleb128 .LVL217-.Ltext0
	.uleb128 .LVL218-.Ltext0
	.uleb128 0x3
	.byte	0x74
	.sleb128 -1
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL218-.Ltext0
	.uleb128 .LVL227-.Ltext0
	.uleb128 0x1
	.byte	0x5d
	.byte	0x4
	.uleb128 .LVL230-.Ltext0
	.uleb128 .LVL231-.Ltext0
	.uleb128 0x1
	.byte	0x5d
	.byte	0x4
	.uleb128 .LVL231-.Ltext0
	.uleb128 .LVL232-.Ltext0
	.uleb128 0x6
	.byte	0xa3
	.uleb128 0x1
	.byte	0x54
	.byte	0x31
	.byte	0x1c
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL232-.Ltext0
	.uleb128 .LVL235-.Ltext0
	.uleb128 0x1
	.byte	0x5d
	.byte	0x4
	.uleb128 .LVL235-.Ltext0
	.uleb128 .LFE68-.Ltext0
	.uleb128 0x3
	.byte	0x74
	.sleb128 -1
	.byte	0x9f
	.byte	0
.LVUS154:
	.uleb128 .LVU745
	.uleb128 .LVU758
	.uleb128 .LVU758
	.uleb128 .LVU778
	.uleb128 .LVU783
	.uleb128 .LVU787
	.uleb128 .LVU787
	.uleb128 0
.LLST154:
	.byte	0x4
	.uleb128 .LVL217-.Ltext0
	.uleb128 .LVL221-1-.Ltext0
	.uleb128 0x1
	.byte	0x55
	.byte	0x4
	.uleb128 .LVL221-1-.Ltext0
	.uleb128 .LVL227-.Ltext0
	.uleb128 0x1
	.byte	0x56
	.byte	0x4
	.uleb128 .LVL230-.Ltext0
	.uleb128 .LVL232-.Ltext0
	.uleb128 0x4
	.byte	0xa3
	.uleb128 0x1
	.byte	0x55
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL232-.Ltext0
	.uleb128 .LFE68-.Ltext0
	.uleb128 0x1
	.byte	0x55
	.byte	0
.LVUS155:
	.uleb128 .LVU747
	.uleb128 .LVU765
	.uleb128 .LVU765
	.uleb128 .LVU778
	.uleb128 .LVU783
	.uleb128 .LVU789
	.uleb128 .LVU789
	.uleb128 .LVU801
	.uleb128 .LVU801
	.uleb128 0
.LLST155:
	.byte	0x4
	.uleb128 .LVL217-.Ltext0
	.uleb128 .LVL224-.Ltext0
	.uleb128 0x2
	.byte	0x30
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL224-.Ltext0
	.uleb128 .LVL227-.Ltext0
	.uleb128 0x1
	.byte	0x53
	.byte	0x4
	.uleb128 .LVL230-.Ltext0
	.uleb128 .LVL233-.Ltext0
	.uleb128 0x2
	.byte	0x30
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL233-.Ltext0
	.uleb128 .LVL234-.Ltext0
	.uleb128 0x1
	.byte	0x53
	.byte	0x4
	.uleb128 .LVL234-.Ltext0
	.uleb128 .LFE68-.Ltext0
	.uleb128 0x1
	.byte	0x50
	.byte	0
.LVUS156:
	.uleb128 .LVU748
	.uleb128 .LVU778
	.uleb128 .LVU783
	.uleb128 0
.LLST156:
	.byte	0x4
	.uleb128 .LVL217-.Ltext0
	.uleb128 .LVL227-.Ltext0
	.uleb128 0x2
	.byte	0x30
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL230-.Ltext0
	.uleb128 .LFE68-.Ltext0
	.uleb128 0x2
	.byte	0x30
	.byte	0x9f
	.byte	0
.LVUS157:
	.uleb128 .LVU749
	.uleb128 .LVU759
	.uleb128 .LVU759
	.uleb128 .LVU763
	.uleb128 .LVU763
	.uleb128 .LVU778
	.uleb128 .LVU787
	.uleb128 0
.LLST157:
	.byte	0x4
	.uleb128 .LVL217-.Ltext0
	.uleb128 .LVL222-.Ltext0
	.uleb128 0x3
	.byte	0x9
	.byte	0xff
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL222-.Ltext0
	.uleb128 .LVL223-.Ltext0
	.uleb128 0x1
	.byte	0x50
	.byte	0x4
	.uleb128 .LVL223-.Ltext0
	.uleb128 .LVL227-.Ltext0
	.uleb128 0x1
	.byte	0x5c
	.byte	0x4
	.uleb128 .LVL232-.Ltext0
	.uleb128 .LFE68-.Ltext0
	.uleb128 0x3
	.byte	0x9
	.byte	0xff
	.byte	0x9f
	.byte	0
.LVUS142:
	.uleb128 0
	.uleb128 .LVU722
	.uleb128 .LVU722
	.uleb128 .LVU734
	.uleb128 .LVU734
	.uleb128 .LVU736
	.uleb128 .LVU736
	.uleb128 0
.LLST142:
	.byte	0x4
	.uleb128 .LVL204-.Ltext0
	.uleb128 .LVL205-.Ltext0
	.uleb128 0x1
	.byte	0x54
	.byte	0x4
	.uleb128 .LVL205-.Ltext0
	.uleb128 .LVL212-.Ltext0
	.uleb128 0x1
	.byte	0x5c
	.byte	0x4
	.uleb128 .LVL212-.Ltext0
	.uleb128 .LVL214-.Ltext0
	.uleb128 0x4
	.byte	0xa3
	.uleb128 0x1
	.byte	0x54
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL214-.Ltext0
	.uleb128 .LFE67-.Ltext0
	.uleb128 0x1
	.byte	0x5c
	.byte	0
.LVUS143:
	.uleb128 0
	.uleb128 .LVU727
	.uleb128 .LVU727
	.uleb128 .LVU735
	.uleb128 .LVU735
	.uleb128 .LVU736
	.uleb128 .LVU736
	.uleb128 0
.LLST143:
	.byte	0x4
	.uleb128 .LVL204-.Ltext0
	.uleb128 .LVL207-1-.Ltext0
	.uleb128 0x1
	.byte	0x51
	.byte	0x4
	.uleb128 .LVL207-1-.Ltext0
	.uleb128 .LVL213-.Ltext0
	.uleb128 0x1
	.byte	0x5d
	.byte	0x4
	.uleb128 .LVL213-.Ltext0
	.uleb128 .LVL214-.Ltext0
	.uleb128 0x4
	.byte	0xa3
	.uleb128 0x1
	.byte	0x51
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL214-.Ltext0
	.uleb128 .LFE67-.Ltext0
	.uleb128 0x1
	.byte	0x5d
	.byte	0
.LVUS144:
	.uleb128 0
	.uleb128 .LVU727
	.uleb128 .LVU727
	.uleb128 0
.LLST144:
	.byte	0x4
	.uleb128 .LVL204-.Ltext0
	.uleb128 .LVL207-1-.Ltext0
	.uleb128 0x1
	.byte	0x52
	.byte	0x4
	.uleb128 .LVL207-1-.Ltext0
	.uleb128 .LFE67-.Ltext0
	.uleb128 0x4
	.byte	0xa3
	.uleb128 0x1
	.byte	0x52
	.byte	0x9f
	.byte	0
.LVUS123:
	.uleb128 0
	.uleb128 .LVU606
	.uleb128 .LVU606
	.uleb128 .LVU674
	.uleb128 .LVU674
	.uleb128 .LVU676
	.uleb128 .LVU676
	.uleb128 0
.LLST123:
	.byte	0x4
	.uleb128 .LVL167-.Ltext0
	.uleb128 .LVL168-.Ltext0
	.uleb128 0x1
	.byte	0x55
	.byte	0x4
	.uleb128 .LVL168-.Ltext0
	.uleb128 .LVL189-.Ltext0
	.uleb128 0x1
	.byte	0x53
	.byte	0x4
	.uleb128 .LVL189-.Ltext0
	.uleb128 .LVL191-.Ltext0
	.uleb128 0x4
	.byte	0xa3
	.uleb128 0x1
	.byte	0x55
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL191-.Ltext0
	.uleb128 .LFE66-.Ltext0
	.uleb128 0x1
	.byte	0x53
	.byte	0
.LVUS124:
	.uleb128 0
	.uleb128 .LVU606
	.uleb128 .LVU606
	.uleb128 0
.LLST124:
	.byte	0x4
	.uleb128 .LVL167-.Ltext0
	.uleb128 .LVL168-.Ltext0
	.uleb128 0x1
	.byte	0x54
	.byte	0x4
	.uleb128 .LVL168-.Ltext0
	.uleb128 .LFE66-.Ltext0
	.uleb128 0x4
	.byte	0xa3
	.uleb128 0x1
	.byte	0x54
	.byte	0x9f
	.byte	0
.LVUS125:
	.uleb128 0
	.uleb128 .LVU606
	.uleb128 .LVU606
	.uleb128 .LVU675
	.uleb128 .LVU675
	.uleb128 .LVU676
	.uleb128 .LVU676
	.uleb128 0
.LLST125:
	.byte	0x4
	.uleb128 .LVL167-.Ltext0
	.uleb128 .LVL168-.Ltext0
	.uleb128 0x1
	.byte	0x51
	.byte	0x4
	.uleb128 .LVL168-.Ltext0
	.uleb128 .LVL190-.Ltext0
	.uleb128 0x1
	.byte	0x5d
	.byte	0x4
	.uleb128 .LVL190-.Ltext0
	.uleb128 .LVL191-.Ltext0
	.uleb128 0x4
	.byte	0xa3
	.uleb128 0x1
	.byte	0x51
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL191-.Ltext0
	.uleb128 .LFE66-.Ltext0
	.uleb128 0x1
	.byte	0x5d
	.byte	0
.LVUS126:
	.uleb128 .LVU646
	.uleb128 .LVU655
.LLST126:
	.byte	0x4
	.uleb128 .LVL182-.Ltext0
	.uleb128 .LVL184-1-.Ltext0
	.uleb128 0x1
	.byte	0x61
	.byte	0
.LVUS127:
	.uleb128 .LVU653
	.uleb128 .LVU655
	.uleb128 .LVU655
	.uleb128 .LVU661
	.uleb128 .LVU676
	.uleb128 .LVU680
	.uleb128 .LVU704
	.uleb128 .LVU717
.LLST127:
	.byte	0x4
	.uleb128 .LVL183-.Ltext0
	.uleb128 .LVL184-1-.Ltext0
	.uleb128 0x1
	.byte	0x62
	.byte	0x4
	.uleb128 .LVL184-1-.Ltext0
	.uleb128 .LVL185-.Ltext0
	.uleb128 0x2
	.byte	0x76
	.sleb128 -56
	.byte	0x4
	.uleb128 .LVL191-.Ltext0
	.uleb128 .LVL192-.Ltext0
	.uleb128 0x2
	.byte	0x76
	.sleb128 -56
	.byte	0x4
	.uleb128 .LVL199-.Ltext0
	.uleb128 .LVL203-.Ltext0
	.uleb128 0x2
	.byte	0x76
	.sleb128 -56
	.byte	0
.LVUS128:
	.uleb128 .LVU656
	.uleb128 .LVU661
	.uleb128 .LVU661
	.uleb128 .LVU662
	.uleb128 .LVU662
	.uleb128 .LVU673
	.uleb128 .LVU676
	.uleb128 .LVU681
	.uleb128 .LVU681
	.uleb128 .LVU688
	.uleb128 .LVU704
	.uleb128 .LVU709
	.uleb128 .LVU709
	.uleb128 .LVU711
	.uleb128 .LVU711
	.uleb128 .LVU715
	.uleb128 .LVU715
	.uleb128 .LVU717
.LLST128:
	.byte	0x4
	.uleb128 .LVL184-.Ltext0
	.uleb128 .LVL185-.Ltext0
	.uleb128 0x15
	.byte	0x76
	.sleb128 -56
	.byte	0xa6
	.byte	0x8
	.uleb128 0x2e
	.byte	0xa4
	.uleb128 0x2e
	.byte	0x8
	.long	0
	.long	0
	.byte	0x2b
	.byte	0x8
	.byte	0xff
	.byte	0x1a
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL185-.Ltext0
	.uleb128 .LVL186-1-.Ltext0
	.uleb128 0x13
	.byte	0xa5
	.uleb128 0x12
	.uleb128 0x2e
	.byte	0xa4
	.uleb128 0x2e
	.byte	0x8
	.long	0
	.long	0
	.byte	0x2b
	.byte	0x8
	.byte	0xff
	.byte	0x1a
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL186-1-.Ltext0
	.uleb128 .LVL188-.Ltext0
	.uleb128 0x15
	.byte	0x76
	.sleb128 -56
	.byte	0xa6
	.byte	0x8
	.uleb128 0x2e
	.byte	0xa4
	.uleb128 0x2e
	.byte	0x8
	.long	0
	.long	0
	.byte	0x2b
	.byte	0x8
	.byte	0xff
	.byte	0x1a
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL191-.Ltext0
	.uleb128 .LVL193-1-.Ltext0
	.uleb128 0x13
	.byte	0xa5
	.uleb128 0x12
	.uleb128 0x2e
	.byte	0xa4
	.uleb128 0x2e
	.byte	0x8
	.long	0
	.long	0
	.byte	0x2b
	.byte	0x8
	.byte	0xff
	.byte	0x1a
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL193-1-.Ltext0
	.uleb128 .LVL194-.Ltext0
	.uleb128 0x15
	.byte	0x76
	.sleb128 -56
	.byte	0xa6
	.byte	0x8
	.uleb128 0x2e
	.byte	0xa4
	.uleb128 0x2e
	.byte	0x8
	.long	0
	.long	0
	.byte	0x2b
	.byte	0x8
	.byte	0xff
	.byte	0x1a
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL199-.Ltext0
	.uleb128 .LVL200-.Ltext0
	.uleb128 0x13
	.byte	0xa5
	.uleb128 0x12
	.uleb128 0x2e
	.byte	0xa4
	.uleb128 0x2e
	.byte	0x8
	.long	0
	.long	0
	.byte	0x2b
	.byte	0x8
	.byte	0xff
	.byte	0x1a
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL200-.Ltext0
	.uleb128 .LVL201-.Ltext0
	.uleb128 0x15
	.byte	0x76
	.sleb128 -56
	.byte	0xa6
	.byte	0x8
	.uleb128 0x2e
	.byte	0xa4
	.uleb128 0x2e
	.byte	0x8
	.long	0
	.long	0
	.byte	0x2b
	.byte	0x8
	.byte	0xff
	.byte	0x1a
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL201-.Ltext0
	.uleb128 .LVL202-.Ltext0
	.uleb128 0x13
	.byte	0xa5
	.uleb128 0x12
	.uleb128 0x2e
	.byte	0xa4
	.uleb128 0x2e
	.byte	0x8
	.long	0
	.long	0
	.byte	0x2b
	.byte	0x8
	.byte	0xff
	.byte	0x1a
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL202-.Ltext0
	.uleb128 .LVL203-.Ltext0
	.uleb128 0x15
	.byte	0x76
	.sleb128 -56
	.byte	0xa6
	.byte	0x8
	.uleb128 0x2e
	.byte	0xa4
	.uleb128 0x2e
	.byte	0x8
	.long	0
	.long	0
	.byte	0x2b
	.byte	0x8
	.byte	0xff
	.byte	0x1a
	.byte	0x9f
	.byte	0
.LVUS129:
	.uleb128 .LVU663
	.uleb128 .LVU669
	.uleb128 .LVU669
	.uleb128 .LVU676
	.uleb128 .LVU682
	.uleb128 .LVU688
	.uleb128 .LVU688
	.uleb128 .LVU700
	.uleb128 .LVU706
	.uleb128 .LVU711
	.uleb128 .LVU712
	.uleb128 .LVU717
.LLST129:
	.byte	0x4
	.uleb128 .LVL186-.Ltext0
	.uleb128 .LVL187-.Ltext0
	.uleb128 0x3
	.byte	0x7d
	.sleb128 16
	.byte	0x6
	.byte	0x4
	.uleb128 .LVL187-.Ltext0
	.uleb128 .LVL191-.Ltext0
	.uleb128 0x1
	.byte	0x61
	.byte	0x4
	.uleb128 .LVL193-.Ltext0
	.uleb128 .LVL194-.Ltext0
	.uleb128 0x3
	.byte	0x7d
	.sleb128 16
	.byte	0x6
	.byte	0x4
	.uleb128 .LVL194-.Ltext0
	.uleb128 .LVL196-.Ltext0
	.uleb128 0x1
	.byte	0x61
	.byte	0x4
	.uleb128 .LVL199-.Ltext0
	.uleb128 .LVL201-.Ltext0
	.uleb128 0x3
	.byte	0x7d
	.sleb128 16
	.byte	0x6
	.byte	0x4
	.uleb128 .LVL201-.Ltext0
	.uleb128 .LVL203-.Ltext0
	.uleb128 0x3
	.byte	0x7d
	.sleb128 16
	.byte	0x6
	.byte	0
.LVUS130:
	.uleb128 .LVU664
	.uleb128 .LVU669
	.uleb128 .LVU683
	.uleb128 .LVU688
	.uleb128 .LVU707
	.uleb128 .LVU709
	.uleb128 .LVU709
	.uleb128 .LVU711
	.uleb128 .LVU713
	.uleb128 .LVU715
	.uleb128 .LVU715
	.uleb128 .LVU717
.LLST130:
	.byte	0x4
	.uleb128 .LVL186-.Ltext0
	.uleb128 .LVL187-.Ltext0
	.uleb128 0x12
	.byte	0x76
	.sleb128 -56
	.byte	0xa6
	.byte	0x8
	.uleb128 0x2e
	.byte	0x12
	.byte	0x1e
	.byte	0x7d
	.sleb128 16
	.byte	0x6
	.byte	0xa6
	.byte	0x8
	.uleb128 0x2e
	.byte	0x2d
	.byte	0x8
	.byte	0xff
	.byte	0x1a
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL193-.Ltext0
	.uleb128 .LVL194-.Ltext0
	.uleb128 0x12
	.byte	0x76
	.sleb128 -56
	.byte	0xa6
	.byte	0x8
	.uleb128 0x2e
	.byte	0x12
	.byte	0x1e
	.byte	0x7d
	.sleb128 16
	.byte	0x6
	.byte	0xa6
	.byte	0x8
	.uleb128 0x2e
	.byte	0x2d
	.byte	0x8
	.byte	0xff
	.byte	0x1a
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL199-.Ltext0
	.uleb128 .LVL200-.Ltext0
	.uleb128 0x10
	.byte	0xa5
	.uleb128 0x12
	.uleb128 0x2e
	.byte	0x12
	.byte	0x1e
	.byte	0x7d
	.sleb128 16
	.byte	0x6
	.byte	0xa6
	.byte	0x8
	.uleb128 0x2e
	.byte	0x2d
	.byte	0x8
	.byte	0xff
	.byte	0x1a
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL200-.Ltext0
	.uleb128 .LVL201-.Ltext0
	.uleb128 0x12
	.byte	0x76
	.sleb128 -56
	.byte	0xa6
	.byte	0x8
	.uleb128 0x2e
	.byte	0x12
	.byte	0x1e
	.byte	0x7d
	.sleb128 16
	.byte	0x6
	.byte	0xa6
	.byte	0x8
	.uleb128 0x2e
	.byte	0x2d
	.byte	0x8
	.byte	0xff
	.byte	0x1a
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL201-.Ltext0
	.uleb128 .LVL202-.Ltext0
	.uleb128 0x10
	.byte	0xa5
	.uleb128 0x12
	.uleb128 0x2e
	.byte	0x12
	.byte	0x1e
	.byte	0x7d
	.sleb128 16
	.byte	0x6
	.byte	0xa6
	.byte	0x8
	.uleb128 0x2e
	.byte	0x2d
	.byte	0x8
	.byte	0xff
	.byte	0x1a
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL202-.Ltext0
	.uleb128 .LVL203-.Ltext0
	.uleb128 0x12
	.byte	0x76
	.sleb128 -56
	.byte	0xa6
	.byte	0x8
	.uleb128 0x2e
	.byte	0x12
	.byte	0x1e
	.byte	0x7d
	.sleb128 16
	.byte	0x6
	.byte	0xa6
	.byte	0x8
	.uleb128 0x2e
	.byte	0x2d
	.byte	0x8
	.byte	0xff
	.byte	0x1a
	.byte	0x9f
	.byte	0
.LVUS132:
	.uleb128 .LVU609
	.uleb128 .LVU646
	.uleb128 .LVU695
	.uleb128 .LVU704
	.uleb128 .LVU717
	.uleb128 0
.LLST132:
	.byte	0x4
	.uleb128 .LVL169-.Ltext0
	.uleb128 .LVL182-.Ltext0
	.uleb128 0x1
	.byte	0x51
	.byte	0x4
	.uleb128 .LVL195-.Ltext0
	.uleb128 .LVL199-.Ltext0
	.uleb128 0x1
	.byte	0x51
	.byte	0x4
	.uleb128 .LVL203-.Ltext0
	.uleb128 .LFE66-.Ltext0
	.uleb128 0x1
	.byte	0x51
	.byte	0
.LVUS133:
	.uleb128 .LVU611
	.uleb128 .LVU616
	.uleb128 .LVU616
	.uleb128 .LVU620
	.uleb128 .LVU622
	.uleb128 .LVU626
	.uleb128 .LVU627
	.uleb128 .LVU631
	.uleb128 .LVU632
	.uleb128 .LVU644
	.uleb128 .LVU697
	.uleb128 .LVU702
	.uleb128 .LVU702
	.uleb128 .LVU704
	.uleb128 .LVU717
	.uleb128 0
.LLST133:
	.byte	0x4
	.uleb128 .LVL169-.Ltext0
	.uleb128 .LVL170-.Ltext0
	.uleb128 0xa
	.byte	0x9e
	.uleb128 0x8
	.long	0
	.long	0
	.byte	0x4
	.uleb128 .LVL170-.Ltext0
	.uleb128 .LVL172-.Ltext0
	.uleb128 0x1
	.byte	0x61
	.byte	0x4
	.uleb128 .LVL173-.Ltext0
	.uleb128 .LVL174-.Ltext0
	.uleb128 0x1
	.byte	0x61
	.byte	0x4
	.uleb128 .LVL175-.Ltext0
	.uleb128 .LVL177-.Ltext0
	.uleb128 0x1
	.byte	0x61
	.byte	0x4
	.uleb128 .LVL178-.Ltext0
	.uleb128 .LVL182-.Ltext0
	.uleb128 0x1
	.byte	0x61
	.byte	0x4
	.uleb128 .LVL195-.Ltext0
	.uleb128 .LVL197-.Ltext0
	.uleb128 0xa
	.byte	0x9e
	.uleb128 0x8
	.long	0
	.long	0
	.byte	0x4
	.uleb128 .LVL197-.Ltext0
	.uleb128 .LVL199-.Ltext0
	.uleb128 0x1
	.byte	0x61
	.byte	0x4
	.uleb128 .LVL203-.Ltext0
	.uleb128 .LFE66-.Ltext0
	.uleb128 0xa
	.byte	0x9e
	.uleb128 0x8
	.long	0
	.long	0
	.byte	0
.LVUS135:
	.uleb128 .LVU613
	.uleb128 .LVU616
	.uleb128 .LVU699
	.uleb128 .LVU702
	.uleb128 .LVU717
	.uleb128 0
.LLST135:
	.byte	0x4
	.uleb128 .LVL169-.Ltext0
	.uleb128 .LVL170-.Ltext0
	.uleb128 0x2
	.byte	0x30
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL195-.Ltext0
	.uleb128 .LVL197-.Ltext0
	.uleb128 0x2
	.byte	0x30
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL203-.Ltext0
	.uleb128 .LFE66-.Ltext0
	.uleb128 0x2
	.byte	0x30
	.byte	0x9f
	.byte	0
.LVUS137:
	.uleb128 .LVU638
	.uleb128 .LVU640
	.uleb128 .LVU640
	.uleb128 .LVU644
.LLST137:
	.byte	0x4
	.uleb128 .LVL180-.Ltext0
	.uleb128 .LVL181-.Ltext0
	.uleb128 0x1
	.byte	0x62
	.byte	0x4
	.uleb128 .LVL181-.Ltext0
	.uleb128 .LVL182-.Ltext0
	.uleb128 0x16
	.byte	0x70
	.sleb128 0
	.byte	0x33
	.byte	0x24
	.byte	0x73
	.sleb128 0
	.byte	0x22
	.byte	0xa6
	.byte	0x8
	.uleb128 0x2e
	.byte	0x70
	.sleb128 0
	.byte	0x33
	.byte	0x24
	.byte	0x71
	.sleb128 0
	.byte	0x22
	.byte	0xa6
	.byte	0x8
	.uleb128 0x2e
	.byte	0x1c
	.byte	0x9f
	.byte	0
.LVUS139:
	.uleb128 .LVU648
	.uleb128 .LVU653
.LLST139:
	.byte	0x4
	.uleb128 .LVL182-.Ltext0
	.uleb128 .LVL183-.Ltext0
	.uleb128 0x2
	.byte	0x7c
	.sleb128 4
	.byte	0
.LVUS140:
	.uleb128 .LVU648
	.uleb128 .LVU653
.LLST140:
	.byte	0x4
	.uleb128 .LVL182-.Ltext0
	.uleb128 .LVL183-.Ltext0
	.uleb128 0x1
	.byte	0x51
	.byte	0
.LVUS141:
	.uleb128 .LVU648
	.uleb128 .LVU653
.LLST141:
	.byte	0x4
	.uleb128 .LVL182-.Ltext0
	.uleb128 .LVL183-.Ltext0
	.uleb128 0x1
	.byte	0x53
	.byte	0
.LVUS30:
	.uleb128 0
	.uleb128 .LVU177
	.uleb128 .LVU177
	.uleb128 .LVU180
	.uleb128 .LVU180
	.uleb128 .LVU241
	.uleb128 .LVU241
	.uleb128 0
.LLST30:
	.byte	0x4
	.uleb128 .LVL48-.Ltext0
	.uleb128 .LVL50-.Ltext0
	.uleb128 0x1
	.byte	0x55
	.byte	0x4
	.uleb128 .LVL50-.Ltext0
	.uleb128 .LVL51-1-.Ltext0
	.uleb128 0x1
	.byte	0x51
	.byte	0x4
	.uleb128 .LVL51-1-.Ltext0
	.uleb128 .LVL64-.Ltext0
	.uleb128 0x1
	.byte	0x56
	.byte	0x4
	.uleb128 .LVL64-.Ltext0
	.uleb128 .LFE60-.Ltext0
	.uleb128 0x4
	.byte	0xa3
	.uleb128 0x1
	.byte	0x55
	.byte	0x9f
	.byte	0
.LVUS35:
	.uleb128 .LVU193
	.uleb128 .LVU194
	.uleb128 .LVU194
	.uleb128 .LVU200
	.uleb128 .LVU200
	.uleb128 .LVU203
.LLST35:
	.byte	0x4
	.uleb128 .LVL53-.Ltext0
	.uleb128 .LVL54-.Ltext0
	.uleb128 0x2
	.byte	0x30
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL54-.Ltext0
	.uleb128 .LVL55-.Ltext0
	.uleb128 0x1
	.byte	0x53
	.byte	0x4
	.uleb128 .LVL55-.Ltext0
	.uleb128 .LVL56-.Ltext0
	.uleb128 0x3
	.byte	0x73
	.sleb128 -1
	.byte	0x9f
	.byte	0
.LVUS37:
	.uleb128 .LVU195
	.uleb128 .LVU201
.LLST37:
	.byte	0x4
	.uleb128 .LVL54-.Ltext0
	.uleb128 .LVL56-.Ltext0
	.uleb128 0xa
	.byte	0x3
	.quad	.LC9
	.byte	0x9f
	.byte	0
.LVUS32:
	.uleb128 .LVU171
	.uleb128 .LVU180
.LLST32:
	.byte	0x4
	.uleb128 .LVL49-.Ltext0
	.uleb128 .LVL51-.Ltext0
	.uleb128 0xa
	.byte	0x3
	.quad	.LC6
	.byte	0x9f
	.byte	0
.LVUS33:
	.uleb128 .LVU182
	.uleb128 .LVU185
.LLST33:
	.byte	0x4
	.uleb128 .LVL51-.Ltext0
	.uleb128 .LVL52-.Ltext0
	.uleb128 0xa
	.byte	0x3
	.quad	.LC7
	.byte	0x9f
	.byte	0
.LVUS34:
	.uleb128 .LVU187
	.uleb128 .LVU190
.LLST34:
	.byte	0x4
	.uleb128 .LVL52-.Ltext0
	.uleb128 .LVL53-.Ltext0
	.uleb128 0xa
	.byte	0x3
	.quad	.LC8
	.byte	0x9f
	.byte	0
.LVUS38:
	.uleb128 .LVU206
	.uleb128 .LVU209
.LLST38:
	.byte	0x4
	.uleb128 .LVL57-.Ltext0
	.uleb128 .LVL58-.Ltext0
	.uleb128 0x6
	.byte	0xa0
	.long	.Ldebug_info0+6107
	.sleb128 0
	.byte	0
.LVUS39:
	.uleb128 .LVU211
	.uleb128 .LVU214
.LLST39:
	.byte	0x4
	.uleb128 .LVL58-.Ltext0
	.uleb128 .LVL59-.Ltext0
	.uleb128 0xa
	.byte	0x3
	.quad	.LC10
	.byte	0x9f
	.byte	0
.LVUS40:
	.uleb128 .LVU216
	.uleb128 .LVU219
.LLST40:
	.byte	0x4
	.uleb128 .LVL59-.Ltext0
	.uleb128 .LVL60-.Ltext0
	.uleb128 0xa
	.byte	0x3
	.quad	.LC11
	.byte	0x9f
	.byte	0
.LVUS41:
	.uleb128 .LVU221
	.uleb128 .LVU224
.LLST41:
	.byte	0x4
	.uleb128 .LVL60-.Ltext0
	.uleb128 .LVL61-.Ltext0
	.uleb128 0xa
	.byte	0x3
	.quad	.LC12
	.byte	0x9f
	.byte	0
.LVUS42:
	.uleb128 .LVU226
	.uleb128 .LVU229
.LLST42:
	.byte	0x4
	.uleb128 .LVL61-.Ltext0
	.uleb128 .LVL62-.Ltext0
	.uleb128 0xa
	.byte	0x3
	.quad	.LC13
	.byte	0x9f
	.byte	0
.LVUS43:
	.uleb128 .LVU231
	.uleb128 .LVU234
.LLST43:
	.byte	0x4
	.uleb128 .LVL62-.Ltext0
	.uleb128 .LVL63-.Ltext0
	.uleb128 0xa
	.byte	0x3
	.quad	.LC14
	.byte	0x9f
	.byte	0
.LVUS26:
	.uleb128 0
	.uleb128 .LVU129
	.uleb128 .LVU129
	.uleb128 .LVU157
	.uleb128 .LVU157
	.uleb128 0
.LLST26:
	.byte	0x4
	.uleb128 .LVL38-.Ltext0
	.uleb128 .LVL40-.Ltext0
	.uleb128 0x1
	.byte	0x51
	.byte	0x4
	.uleb128 .LVL40-.Ltext0
	.uleb128 .LVL45-.Ltext0
	.uleb128 0x1
	.byte	0x52
	.byte	0x4
	.uleb128 .LVL45-.Ltext0
	.uleb128 .LFE58-.Ltext0
	.uleb128 0x1
	.byte	0x51
	.byte	0
.LVUS28:
	.uleb128 .LVU126
	.uleb128 .LVU130
	.uleb128 .LVU157
	.uleb128 0
.LLST28:
	.byte	0x4
	.uleb128 .LVL39-.Ltext0
	.uleb128 .LVL41-.Ltext0
	.uleb128 0x2
	.byte	0x30
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL45-.Ltext0
	.uleb128 .LFE58-.Ltext0
	.uleb128 0x2
	.byte	0x30
	.byte	0x9f
	.byte	0
.LVUS23:
	.uleb128 0
	.uleb128 .LVU103
	.uleb128 .LVU103
	.uleb128 0
.LLST23:
	.byte	0x4
	.uleb128 .LVL31-.Ltext0
	.uleb128 .LVL33-.Ltext0
	.uleb128 0x1
	.byte	0x55
	.byte	0x4
	.uleb128 .LVL33-.Ltext0
	.uleb128 .LFE57-.Ltext0
	.uleb128 0x4
	.byte	0xa3
	.uleb128 0x1
	.byte	0x55
	.byte	0x9f
	.byte	0
.LVUS24:
	.uleb128 0
	.uleb128 .LVU103
	.uleb128 .LVU103
	.uleb128 0
.LLST24:
	.byte	0x4
	.uleb128 .LVL31-.Ltext0
	.uleb128 .LVL33-.Ltext0
	.uleb128 0x1
	.byte	0x54
	.byte	0x4
	.uleb128 .LVL33-.Ltext0
	.uleb128 .LFE57-.Ltext0
	.uleb128 0x4
	.byte	0xa3
	.uleb128 0x1
	.byte	0x54
	.byte	0x9f
	.byte	0
.LVUS25:
	.uleb128 .LVU99
	.uleb128 .LVU103
	.uleb128 .LVU103
	.uleb128 .LVU106
	.uleb128 .LVU106
	.uleb128 .LVU107
	.uleb128 .LVU107
	.uleb128 .LVU119
	.uleb128 .LVU119
	.uleb128 .LVU120
.LLST25:
	.byte	0x4
	.uleb128 .LVL32-.Ltext0
	.uleb128 .LVL33-.Ltext0
	.uleb128 0x2
	.byte	0x30
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL33-.Ltext0
	.uleb128 .LVL34-.Ltext0
	.uleb128 0x1
	.byte	0x50
	.byte	0x4
	.uleb128 .LVL34-.Ltext0
	.uleb128 .LVL35-.Ltext0
	.uleb128 0x2
	.byte	0x75
	.sleb128 16
	.byte	0x4
	.uleb128 .LVL35-.Ltext0
	.uleb128 .LVL36-.Ltext0
	.uleb128 0x2
	.byte	0x75
	.sleb128 -32
	.byte	0x4
	.uleb128 .LVL36-.Ltext0
	.uleb128 .LVL37-.Ltext0
	.uleb128 0x1
	.byte	0x50
	.byte	0
.LVUS15:
	.uleb128 .LVU68
	.uleb128 .LVU71
.LLST15:
	.byte	0x4
	.uleb128 .LVL25-.Ltext0
	.uleb128 .LVL26-.Ltext0
	.uleb128 0x2
	.byte	0x40
	.byte	0x9f
	.byte	0
.LVUS16:
	.uleb128 .LVU68
	.uleb128 .LVU71
.LLST16:
	.byte	0x4
	.uleb128 .LVL25-.Ltext0
	.uleb128 .LVL26-.Ltext0
	.uleb128 0x1
	.byte	0x55
	.byte	0
.LVUS17:
	.uleb128 .LVU73
	.uleb128 .LVU76
.LLST17:
	.byte	0x4
	.uleb128 .LVL26-.Ltext0
	.uleb128 .LVL27-.Ltext0
	.uleb128 0x2
	.byte	0x40
	.byte	0x9f
	.byte	0
.LVUS18:
	.uleb128 .LVU73
	.uleb128 .LVU76
.LLST18:
	.byte	0x4
	.uleb128 .LVL26-.Ltext0
	.uleb128 .LVL27-.Ltext0
	.uleb128 0x1
	.byte	0x54
	.byte	0
.LVUS19:
	.uleb128 .LVU73
	.uleb128 .LVU76
.LLST19:
	.byte	0x4
	.uleb128 .LVL26-.Ltext0
	.uleb128 .LVL27-.Ltext0
	.uleb128 0x1
	.byte	0x55
	.byte	0
.LVUS20:
	.uleb128 .LVU78
	.uleb128 .LVU81
.LLST20:
	.byte	0x4
	.uleb128 .LVL27-.Ltext0
	.uleb128 .LVL28-.Ltext0
	.uleb128 0x2
	.byte	0x40
	.byte	0x9f
	.byte	0
.LVUS21:
	.uleb128 .LVU78
	.uleb128 .LVU81
.LLST21:
	.byte	0x4
	.uleb128 .LVL27-.Ltext0
	.uleb128 .LVL28-.Ltext0
	.uleb128 0x1
	.byte	0x54
	.byte	0
.LVUS0:
	.uleb128 .LVU4
	.uleb128 .LVU7
.LLST0:
	.byte	0x4
	.uleb128 .LVL1-.Ltext0
	.uleb128 .LVL2-.Ltext0
	.uleb128 0x2
	.byte	0x38
	.byte	0x9f
	.byte	0
.LVUS1:
	.uleb128 .LVU4
	.uleb128 .LVU7
.LLST1:
	.byte	0x4
	.uleb128 .LVL1-.Ltext0
	.uleb128 .LVL2-.Ltext0
	.uleb128 0x1
	.byte	0x55
	.byte	0
.LVUS2:
	.uleb128 .LVU4
	.uleb128 .LVU7
.LLST2:
	.byte	0x4
	.uleb128 .LVL1-.Ltext0
	.uleb128 .LVL2-.Ltext0
	.uleb128 0x6
	.byte	0xa0
	.long	.Ldebug_info0+3649
	.sleb128 0
	.byte	0
.LVUS3:
	.uleb128 .LVU9
	.uleb128 .LVU12
.LLST3:
	.byte	0x4
	.uleb128 .LVL2-.Ltext0
	.uleb128 .LVL3-.Ltext0
	.uleb128 0x2
	.byte	0x38
	.byte	0x9f
	.byte	0
.LVUS4:
	.uleb128 .LVU9
	.uleb128 .LVU12
.LLST4:
	.byte	0x4
	.uleb128 .LVL2-.Ltext0
	.uleb128 .LVL3-.Ltext0
	.uleb128 0x1
	.byte	0x54
	.byte	0
.LVUS5:
	.uleb128 .LVU9
	.uleb128 .LVU12
.LLST5:
	.byte	0x4
	.uleb128 .LVL2-.Ltext0
	.uleb128 .LVL3-.Ltext0
	.uleb128 0x1
	.byte	0x55
	.byte	0
.LVUS6:
	.uleb128 .LVU14
	.uleb128 .LVU17
.LLST6:
	.byte	0x4
	.uleb128 .LVL3-.Ltext0
	.uleb128 .LVL4-.Ltext0
	.uleb128 0x2
	.byte	0x38
	.byte	0x9f
	.byte	0
.LVUS7:
	.uleb128 .LVU14
	.uleb128 .LVU17
.LLST7:
	.byte	0x4
	.uleb128 .LVL3-.Ltext0
	.uleb128 .LVL4-.Ltext0
	.uleb128 0x6
	.byte	0xa0
	.long	.Ldebug_info0+3649
	.sleb128 0
	.byte	0
.LVUS8:
	.uleb128 .LVU14
	.uleb128 .LVU17
.LLST8:
	.byte	0x4
	.uleb128 .LVL3-.Ltext0
	.uleb128 .LVL4-.Ltext0
	.uleb128 0x1
	.byte	0x54
	.byte	0
.LVUS9:
	.uleb128 0
	.uleb128 .LVU28
	.uleb128 .LVU28
	.uleb128 0
.LLST9:
	.byte	0x4
	.uleb128 .LVL5-.Ltext0
	.uleb128 .LVL7-.Ltext0
	.uleb128 0x1
	.byte	0x55
	.byte	0x4
	.uleb128 .LVL7-.Ltext0
	.uleb128 .LFE54-.Ltext0
	.uleb128 0x1
	.byte	0x51
	.byte	0
.LVUS10:
	.uleb128 .LVU23
	.uleb128 .LVU30
	.uleb128 .LVU30
	.uleb128 .LVU34
	.uleb128 .LVU36
	.uleb128 .LVU39
	.uleb128 .LVU40
	.uleb128 .LVU44
	.uleb128 .LVU45
	.uleb128 .LVU59
	.uleb128 .LVU59
	.uleb128 0
.LLST10:
	.byte	0x4
	.uleb128 .LVL6-.Ltext0
	.uleb128 .LVL8-.Ltext0
	.uleb128 0xa
	.byte	0x9e
	.uleb128 0x8
	.long	0
	.long	0
	.byte	0x4
	.uleb128 .LVL8-.Ltext0
	.uleb128 .LVL10-.Ltext0
	.uleb128 0x1
	.byte	0x61
	.byte	0x4
	.uleb128 .LVL11-.Ltext0
	.uleb128 .LVL12-.Ltext0
	.uleb128 0x1
	.byte	0x61
	.byte	0x4
	.uleb128 .LVL13-.Ltext0
	.uleb128 .LVL15-.Ltext0
	.uleb128 0x1
	.byte	0x61
	.byte	0x4
	.uleb128 .LVL16-.Ltext0
	.uleb128 .LVL23-.Ltext0
	.uleb128 0x1
	.byte	0x61
	.byte	0x4
	.uleb128 .LVL23-.Ltext0
	.uleb128 .LFE54-.Ltext0
	.uleb128 0xa
	.byte	0x9e
	.uleb128 0x8
	.long	0
	.long	0
	.byte	0
.LVUS12:
	.uleb128 .LVU25
	.uleb128 .LVU30
	.uleb128 .LVU59
	.uleb128 0
.LLST12:
	.byte	0x4
	.uleb128 .LVL6-.Ltext0
	.uleb128 .LVL8-.Ltext0
	.uleb128 0x2
	.byte	0x30
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL23-.Ltext0
	.uleb128 .LFE54-.Ltext0
	.uleb128 0x2
	.byte	0x30
	.byte	0x9f
	.byte	0
.LVUS14:
	.uleb128 .LVU50
	.uleb128 .LVU52
	.uleb128 .LVU52
	.uleb128 .LVU56
.LLST14:
	.byte	0x4
	.uleb128 .LVL18-.Ltext0
	.uleb128 .LVL19-.Ltext0
	.uleb128 0x1
	.byte	0x62
	.byte	0x4
	.uleb128 .LVL19-.Ltext0
	.uleb128 .LVL21-.Ltext0
	.uleb128 0x16
	.byte	0x70
	.sleb128 0
	.byte	0x33
	.byte	0x24
	.byte	0x71
	.sleb128 0
	.byte	0x22
	.byte	0xa6
	.byte	0x8
	.uleb128 0x2e
	.byte	0x70
	.sleb128 0
	.byte	0x33
	.byte	0x24
	.byte	0x74
	.sleb128 0
	.byte	0x22
	.byte	0xa6
	.byte	0x8
	.uleb128 0x2e
	.byte	0x1c
	.byte	0x9f
	.byte	0
.LVUS22:
	.uleb128 .LVU89
	.uleb128 0
.LLST22:
	.byte	0x4
	.uleb128 .LVL30-.Ltext0
	.uleb128 .LFE56-.Ltext0
	.uleb128 0x1
	.byte	0x50
	.byte	0
.LVUS29:
	.uleb128 .LVU163
	.uleb128 0
.LLST29:
	.byte	0x4
	.uleb128 .LVL47-.Ltext0
	.uleb128 .LFE59-.Ltext0
	.uleb128 0x18
	.byte	0x71
	.sleb128 0
	.byte	0x33
	.byte	0x24
	.byte	0x75
	.sleb128 8
	.byte	0x6
	.byte	0x22
	.byte	0xa6
	.byte	0x8
	.uleb128 0x2e
	.byte	0x71
	.sleb128 0
	.byte	0x33
	.byte	0x24
	.byte	0x74
	.sleb128 8
	.byte	0x6
	.byte	0x22
	.byte	0xa6
	.byte	0x8
	.uleb128 0x2e
	.byte	0x1c
	.byte	0x9f
	.byte	0
.LVUS45:
	.uleb128 0
	.uleb128 .LVU248
	.uleb128 .LVU248
	.uleb128 .LVU263
	.uleb128 .LVU263
	.uleb128 0
.LLST45:
	.byte	0x4
	.uleb128 .LVL66-.Ltext0
	.uleb128 .LVL67-.Ltext0
	.uleb128 0x1
	.byte	0x54
	.byte	0x4
	.uleb128 .LVL67-.Ltext0
	.uleb128 .LVL72-.Ltext0
	.uleb128 0x1
	.byte	0x50
	.byte	0x4
	.uleb128 .LVL72-.Ltext0
	.uleb128 .LFE61-.Ltext0
	.uleb128 0x4
	.byte	0xa3
	.uleb128 0x1
	.byte	0x54
	.byte	0x9f
	.byte	0
.LVUS46:
	.uleb128 0
	.uleb128 .LVU258
	.uleb128 .LVU258
	.uleb128 0
.LLST46:
	.byte	0x4
	.uleb128 .LVL66-.Ltext0
	.uleb128 .LVL70-.Ltext0
	.uleb128 0x1
	.byte	0x51
	.byte	0x4
	.uleb128 .LVL70-.Ltext0
	.uleb128 .LFE61-.Ltext0
	.uleb128 0x4
	.byte	0xa3
	.uleb128 0x1
	.byte	0x51
	.byte	0x9f
	.byte	0
.LVUS47:
	.uleb128 0
	.uleb128 .LVU262
	.uleb128 .LVU262
	.uleb128 0
.LLST47:
	.byte	0x4
	.uleb128 .LVL66-.Ltext0
	.uleb128 .LVL71-.Ltext0
	.uleb128 0x1
	.byte	0x52
	.byte	0x4
	.uleb128 .LVL71-.Ltext0
	.uleb128 .LFE61-.Ltext0
	.uleb128 0x4
	.byte	0xa3
	.uleb128 0x1
	.byte	0x52
	.byte	0x9f
	.byte	0
.LVUS48:
	.uleb128 .LVU251
	.uleb128 .LVU256
	.uleb128 .LVU256
	.uleb128 .LVU263
	.uleb128 .LVU290
	.uleb128 0
.LLST48:
	.byte	0x4
	.uleb128 .LVL68-.Ltext0
	.uleb128 .LVL69-.Ltext0
	.uleb128 0x1
	.byte	0x58
	.byte	0x4
	.uleb128 .LVL69-.Ltext0
	.uleb128 .LVL72-.Ltext0
	.uleb128 0x2
	.byte	0x7a
	.sleb128 0
	.byte	0x4
	.uleb128 .LVL78-.Ltext0
	.uleb128 .LFE61-.Ltext0
	.uleb128 0x1
	.byte	0x58
	.byte	0
.LVUS49:
	.uleb128 .LVU252
	.uleb128 .LVU273
	.uleb128 .LVU273
	.uleb128 .LVU274
	.uleb128 .LVU274
	.uleb128 .LVU285
	.uleb128 .LVU285
	.uleb128 0
.LLST49:
	.byte	0x4
	.uleb128 .LVL68-.Ltext0
	.uleb128 .LVL74-.Ltext0
	.uleb128 0x1
	.byte	0x54
	.byte	0x4
	.uleb128 .LVL74-.Ltext0
	.uleb128 .LVL74-.Ltext0
	.uleb128 0x3
	.byte	0x74
	.sleb128 -1
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL74-.Ltext0
	.uleb128 .LVL77-.Ltext0
	.uleb128 0x3
	.byte	0x74
	.sleb128 1
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL77-.Ltext0
	.uleb128 .LFE61-.Ltext0
	.uleb128 0x1
	.byte	0x54
	.byte	0
.LVUS50:
	.uleb128 .LVU254
	.uleb128 .LVU263
.LLST50:
	.byte	0x4
	.uleb128 .LVL68-.Ltext0
	.uleb128 .LVL72-.Ltext0
	.uleb128 0x1
	.byte	0x50
	.byte	0
.LVUS52:
	.uleb128 .LVU266
	.uleb128 .LVU269
.LLST52:
	.byte	0x4
	.uleb128 .LVL73-.Ltext0
	.uleb128 .LVL73-.Ltext0
	.uleb128 0x4
	.byte	0xa3
	.uleb128 0x1
	.byte	0x52
	.byte	0x9f
	.byte	0
.LVUS53:
	.uleb128 .LVU266
	.uleb128 .LVU269
.LLST53:
	.byte	0x4
	.uleb128 .LVL73-.Ltext0
	.uleb128 .LVL73-.Ltext0
	.uleb128 0x1
	.byte	0x50
	.byte	0
.LVUS54:
	.uleb128 .LVU268
	.uleb128 .LVU269
.LLST54:
	.byte	0x4
	.uleb128 .LVL73-.Ltext0
	.uleb128 .LVL73-.Ltext0
	.uleb128 0xe
	.byte	0x70
	.sleb128 8
	.byte	0x6
	.byte	0x79
	.sleb128 0
	.byte	0x22
	.byte	0xa6
	.byte	0x8
	.uleb128 0x2e
	.byte	0xa5
	.uleb128 0x13
	.uleb128 0x2e
	.byte	0x1c
	.byte	0x9f
	.byte	0
.LVUS55:
	.uleb128 .LVU277
	.uleb128 .LVU285
.LLST55:
	.byte	0x4
	.uleb128 .LVL75-.Ltext0
	.uleb128 .LVL77-.Ltext0
	.uleb128 0x1
	.byte	0x52
	.byte	0
.LVUS56:
	.uleb128 .LVU277
	.uleb128 .LVU285
.LLST56:
	.byte	0x4
	.uleb128 .LVL75-.Ltext0
	.uleb128 .LVL77-.Ltext0
	.uleb128 0x1
	.byte	0x51
	.byte	0
.LVUS57:
	.uleb128 .LVU281
	.uleb128 .LVU285
.LLST57:
	.byte	0x4
	.uleb128 .LVL76-.Ltext0
	.uleb128 .LVL77-.Ltext0
	.uleb128 0x1
	.byte	0x5b
	.byte	0
.LVUS59:
	.uleb128 .LVU293
	.uleb128 .LVU301
.LLST59:
	.byte	0x4
	.uleb128 .LVL79-.Ltext0
	.uleb128 .LVL81-.Ltext0
	.uleb128 0x1
	.byte	0x5a
	.byte	0
.LVUS60:
	.uleb128 .LVU293
	.uleb128 .LVU301
.LLST60:
	.byte	0x4
	.uleb128 .LVL79-.Ltext0
	.uleb128 .LVL81-.Ltext0
	.uleb128 0x1
	.byte	0x51
	.byte	0
.LVUS61:
	.uleb128 .LVU297
	.uleb128 .LVU301
.LLST61:
	.byte	0x4
	.uleb128 .LVL80-.Ltext0
	.uleb128 .LVL81-.Ltext0
	.uleb128 0x1
	.byte	0x52
	.byte	0
.LVUS62:
	.uleb128 0
	.uleb128 .LVU316
	.uleb128 .LVU316
	.uleb128 .LVU405
	.uleb128 .LVU405
	.uleb128 .LVU415
	.uleb128 .LVU415
	.uleb128 .LVU421
	.uleb128 .LVU421
	.uleb128 .LVU434
	.uleb128 .LVU434
	.uleb128 0
.LLST62:
	.byte	0x4
	.uleb128 .LVL82-.Ltext0
	.uleb128 .LVL84-.Ltext0
	.uleb128 0x1
	.byte	0x55
	.byte	0x4
	.uleb128 .LVL84-.Ltext0
	.uleb128 .LVL106-.Ltext0
	.uleb128 0x1
	.byte	0x59
	.byte	0x4
	.uleb128 .LVL106-.Ltext0
	.uleb128 .LVL110-.Ltext0
	.uleb128 0x1
	.byte	0x55
	.byte	0x4
	.uleb128 .LVL110-.Ltext0
	.uleb128 .LVL113-.Ltext0
	.uleb128 0x1
	.byte	0x59
	.byte	0x4
	.uleb128 .LVL113-.Ltext0
	.uleb128 .LVL117-.Ltext0
	.uleb128 0x4
	.byte	0xa3
	.uleb128 0x1
	.byte	0x55
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL117-.Ltext0
	.uleb128 .LFE62-.Ltext0
	.uleb128 0x1
	.byte	0x55
	.byte	0
.LVUS63:
	.uleb128 0
	.uleb128 .LVU326
	.uleb128 .LVU326
	.uleb128 .LVU405
	.uleb128 .LVU405
	.uleb128 .LVU413
	.uleb128 .LVU413
	.uleb128 .LVU433
	.uleb128 .LVU433
	.uleb128 .LVU434
	.uleb128 .LVU434
	.uleb128 0
.LLST63:
	.byte	0x4
	.uleb128 .LVL82-.Ltext0
	.uleb128 .LVL85-.Ltext0
	.uleb128 0x1
	.byte	0x54
	.byte	0x4
	.uleb128 .LVL85-.Ltext0
	.uleb128 .LVL106-.Ltext0
	.uleb128 0x1
	.byte	0x5a
	.byte	0x4
	.uleb128 .LVL106-.Ltext0
	.uleb128 .LVL109-.Ltext0
	.uleb128 0x1
	.byte	0x54
	.byte	0x4
	.uleb128 .LVL109-.Ltext0
	.uleb128 .LVL116-.Ltext0
	.uleb128 0x1
	.byte	0x5a
	.byte	0x4
	.uleb128 .LVL116-.Ltext0
	.uleb128 .LVL117-.Ltext0
	.uleb128 0x4
	.byte	0xa3
	.uleb128 0x1
	.byte	0x54
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL117-.Ltext0
	.uleb128 .LFE62-.Ltext0
	.uleb128 0x1
	.byte	0x54
	.byte	0
.LVUS64:
	.uleb128 0
	.uleb128 .LVU326
	.uleb128 .LVU326
	.uleb128 .LVU400
	.uleb128 .LVU400
	.uleb128 .LVU403
	.uleb128 .LVU405
	.uleb128 .LVU409
	.uleb128 .LVU409
	.uleb128 .LVU420
	.uleb128 .LVU420
	.uleb128 .LVU433
	.uleb128 .LVU433
	.uleb128 .LVU434
	.uleb128 .LVU434
	.uleb128 0
.LLST64:
	.byte	0x4
	.uleb128 .LVL82-.Ltext0
	.uleb128 .LVL85-.Ltext0
	.uleb128 0x1
	.byte	0x51
	.byte	0x4
	.uleb128 .LVL85-.Ltext0
	.uleb128 .LVL104-.Ltext0
	.uleb128 0x1
	.byte	0x5b
	.byte	0x4
	.uleb128 .LVL104-.Ltext0
	.uleb128 .LVL105-.Ltext0
	.uleb128 0x1
	.byte	0x54
	.byte	0x4
	.uleb128 .LVL106-.Ltext0
	.uleb128 .LVL108-.Ltext0
	.uleb128 0x1
	.byte	0x51
	.byte	0x4
	.uleb128 .LVL108-.Ltext0
	.uleb128 .LVL112-.Ltext0
	.uleb128 0x1
	.byte	0x58
	.byte	0x4
	.uleb128 .LVL112-.Ltext0
	.uleb128 .LVL116-.Ltext0
	.uleb128 0x1
	.byte	0x5b
	.byte	0x4
	.uleb128 .LVL116-.Ltext0
	.uleb128 .LVL117-.Ltext0
	.uleb128 0x4
	.byte	0xa3
	.uleb128 0x1
	.byte	0x51
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL117-.Ltext0
	.uleb128 .LFE62-.Ltext0
	.uleb128 0x1
	.byte	0x51
	.byte	0
.LVUS65:
	.uleb128 0
	.uleb128 .LVU326
	.uleb128 .LVU326
	.uleb128 .LVU405
	.uleb128 .LVU405
	.uleb128 .LVU433
	.uleb128 .LVU433
	.uleb128 .LVU434
	.uleb128 .LVU434
	.uleb128 0
.LLST65:
	.byte	0x4
	.uleb128 .LVL82-.Ltext0
	.uleb128 .LVL85-.Ltext0
	.uleb128 0x1
	.byte	0x52
	.byte	0x4
	.uleb128 .LVL85-.Ltext0
	.uleb128 .LVL106-.Ltext0
	.uleb128 0x4
	.byte	0xa3
	.uleb128 0x1
	.byte	0x52
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL106-.Ltext0
	.uleb128 .LVL116-.Ltext0
	.uleb128 0x1
	.byte	0x52
	.byte	0x4
	.uleb128 .LVL116-.Ltext0
	.uleb128 .LVL117-.Ltext0
	.uleb128 0x4
	.byte	0xa3
	.uleb128 0x1
	.byte	0x52
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL117-.Ltext0
	.uleb128 .LFE62-.Ltext0
	.uleb128 0x1
	.byte	0x52
	.byte	0
.LVUS66:
	.uleb128 .LVU307
	.uleb128 .LVU326
	.uleb128 .LVU326
	.uleb128 .LVU405
	.uleb128 .LVU405
	.uleb128 .LVU409
	.uleb128 .LVU409
	.uleb128 .LVU413
	.uleb128 .LVU413
	.uleb128 .LVU420
	.uleb128 .LVU420
	.uleb128 .LVU433
	.uleb128 .LVU433
	.uleb128 .LVU434
	.uleb128 .LVU434
	.uleb128 0
.LLST66:
	.byte	0x4
	.uleb128 .LVL83-.Ltext0
	.uleb128 .LVL85-.Ltext0
	.uleb128 0xd
	.byte	0x71
	.sleb128 0
	.byte	0x74
	.sleb128 0
	.byte	0x1c
	.byte	0x23
	.uleb128 0x1
	.byte	0x32
	.byte	0x1b
	.byte	0x74
	.sleb128 0
	.byte	0x22
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL85-.Ltext0
	.uleb128 .LVL106-.Ltext0
	.uleb128 0x10
	.byte	0xa3
	.uleb128 0x1
	.byte	0x51
	.byte	0xa3
	.uleb128 0x1
	.byte	0x54
	.byte	0x1c
	.byte	0x23
	.uleb128 0x1
	.byte	0x32
	.byte	0x1b
	.byte	0xa3
	.uleb128 0x1
	.byte	0x54
	.byte	0x22
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL106-.Ltext0
	.uleb128 .LVL108-.Ltext0
	.uleb128 0xd
	.byte	0x71
	.sleb128 0
	.byte	0x74
	.sleb128 0
	.byte	0x1c
	.byte	0x23
	.uleb128 0x1
	.byte	0x32
	.byte	0x1b
	.byte	0x74
	.sleb128 0
	.byte	0x22
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL108-.Ltext0
	.uleb128 .LVL109-.Ltext0
	.uleb128 0xd
	.byte	0x78
	.sleb128 0
	.byte	0x74
	.sleb128 0
	.byte	0x1c
	.byte	0x23
	.uleb128 0x1
	.byte	0x32
	.byte	0x1b
	.byte	0x74
	.sleb128 0
	.byte	0x22
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL109-.Ltext0
	.uleb128 .LVL112-.Ltext0
	.uleb128 0xd
	.byte	0x78
	.sleb128 0
	.byte	0x7a
	.sleb128 0
	.byte	0x1c
	.byte	0x23
	.uleb128 0x1
	.byte	0x32
	.byte	0x1b
	.byte	0x7a
	.sleb128 0
	.byte	0x22
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL112-.Ltext0
	.uleb128 .LVL116-.Ltext0
	.uleb128 0xd
	.byte	0x7b
	.sleb128 0
	.byte	0x7a
	.sleb128 0
	.byte	0x1c
	.byte	0x23
	.uleb128 0x1
	.byte	0x32
	.byte	0x1b
	.byte	0x7a
	.sleb128 0
	.byte	0x22
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL116-.Ltext0
	.uleb128 .LVL117-.Ltext0
	.uleb128 0x10
	.byte	0xa3
	.uleb128 0x1
	.byte	0x51
	.byte	0xa3
	.uleb128 0x1
	.byte	0x54
	.byte	0x1c
	.byte	0x23
	.uleb128 0x1
	.byte	0x32
	.byte	0x1b
	.byte	0xa3
	.uleb128 0x1
	.byte	0x54
	.byte	0x22
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL117-.Ltext0
	.uleb128 .LFE62-.Ltext0
	.uleb128 0xd
	.byte	0x71
	.sleb128 0
	.byte	0x74
	.sleb128 0
	.byte	0x1c
	.byte	0x23
	.uleb128 0x1
	.byte	0x32
	.byte	0x1b
	.byte	0x74
	.sleb128 0
	.byte	0x22
	.byte	0x9f
	.byte	0
.LVUS68:
	.uleb128 .LVU380
	.uleb128 .LVU389
	.uleb128 .LVU389
	.uleb128 .LVU399
	.uleb128 .LVU400
	.uleb128 .LVU405
.LLST68:
	.byte	0x4
	.uleb128 .LVL98-.Ltext0
	.uleb128 .LVL100-.Ltext0
	.uleb128 0x3
	.byte	0x74
	.sleb128 1
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL100-.Ltext0
	.uleb128 .LVL103-.Ltext0
	.uleb128 0x1
	.byte	0x5c
	.byte	0x4
	.uleb128 .LVL104-.Ltext0
	.uleb128 .LVL106-.Ltext0
	.uleb128 0x3
	.byte	0x74
	.sleb128 1
	.byte	0x9f
	.byte	0
.LVUS70:
	.uleb128 .LVU327
	.uleb128 .LVU380
	.uleb128 .LVU390
	.uleb128 .LVU400
.LLST70:
	.byte	0x4
	.uleb128 .LVL85-.Ltext0
	.uleb128 .LVL98-.Ltext0
	.uleb128 0x4
	.byte	0xa3
	.uleb128 0x1
	.byte	0x52
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL100-.Ltext0
	.uleb128 .LVL104-.Ltext0
	.uleb128 0x4
	.byte	0xa3
	.uleb128 0x1
	.byte	0x52
	.byte	0x9f
	.byte	0
.LVUS71:
	.uleb128 .LVU327
	.uleb128 .LVU380
	.uleb128 .LVU390
	.uleb128 .LVU400
.LLST71:
	.byte	0x4
	.uleb128 .LVL85-.Ltext0
	.uleb128 .LVL98-.Ltext0
	.uleb128 0x1
	.byte	0x5b
	.byte	0x4
	.uleb128 .LVL100-.Ltext0
	.uleb128 .LVL104-.Ltext0
	.uleb128 0x1
	.byte	0x5b
	.byte	0
.LVUS72:
	.uleb128 .LVU327
	.uleb128 .LVU380
	.uleb128 .LVU390
	.uleb128 .LVU400
.LLST72:
	.byte	0x4
	.uleb128 .LVL85-.Ltext0
	.uleb128 .LVL98-.Ltext0
	.uleb128 0x1
	.byte	0x5a
	.byte	0x4
	.uleb128 .LVL100-.Ltext0
	.uleb128 .LVL104-.Ltext0
	.uleb128 0x1
	.byte	0x5a
	.byte	0
.LVUS73:
	.uleb128 .LVU327
	.uleb128 .LVU380
	.uleb128 .LVU390
	.uleb128 .LVU400
.LLST73:
	.byte	0x4
	.uleb128 .LVL85-.Ltext0
	.uleb128 .LVL98-.Ltext0
	.uleb128 0x1
	.byte	0x59
	.byte	0x4
	.uleb128 .LVL100-.Ltext0
	.uleb128 .LVL104-.Ltext0
	.uleb128 0x1
	.byte	0x59
	.byte	0
.LVUS74:
	.uleb128 .LVU332
	.uleb128 .LVU338
	.uleb128 .LVU364
	.uleb128 .LVU367
	.uleb128 .LVU367
	.uleb128 .LVU380
	.uleb128 .LVU394
	.uleb128 .LVU400
.LLST74:
	.byte	0x4
	.uleb128 .LVL86-.Ltext0
	.uleb128 .LVL87-.Ltext0
	.uleb128 0x1
	.byte	0x50
	.byte	0x4
	.uleb128 .LVL94-.Ltext0
	.uleb128 .LVL95-.Ltext0
	.uleb128 0x2
	.byte	0x73
	.sleb128 0
	.byte	0x4
	.uleb128 .LVL95-.Ltext0
	.uleb128 .LVL98-.Ltext0
	.uleb128 0x1
	.byte	0x50
	.byte	0x4
	.uleb128 .LVL101-.Ltext0
	.uleb128 .LVL104-.Ltext0
	.uleb128 0x1
	.byte	0x50
	.byte	0
.LVUS75:
	.uleb128 .LVU333
	.uleb128 .LVU349
	.uleb128 .LVU349
	.uleb128 .LVU350
	.uleb128 .LVU350
	.uleb128 .LVU361
	.uleb128 .LVU361
	.uleb128 .LVU380
	.uleb128 .LVU396
	.uleb128 .LVU399
	.uleb128 .LVU399
	.uleb128 .LVU400
.LLST75:
	.byte	0x4
	.uleb128 .LVL86-.Ltext0
	.uleb128 .LVL90-.Ltext0
	.uleb128 0x1
	.byte	0x54
	.byte	0x4
	.uleb128 .LVL90-.Ltext0
	.uleb128 .LVL90-.Ltext0
	.uleb128 0x3
	.byte	0x74
	.sleb128 -1
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL90-.Ltext0
	.uleb128 .LVL93-.Ltext0
	.uleb128 0x3
	.byte	0x74
	.sleb128 1
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL93-.Ltext0
	.uleb128 .LVL98-.Ltext0
	.uleb128 0x1
	.byte	0x54
	.byte	0x4
	.uleb128 .LVL102-.Ltext0
	.uleb128 .LVL103-.Ltext0
	.uleb128 0x1
	.byte	0x5c
	.byte	0x4
	.uleb128 .LVL103-.Ltext0
	.uleb128 .LVL104-.Ltext0
	.uleb128 0x1
	.byte	0x54
	.byte	0
.LVUS77:
	.uleb128 .LVU335
	.uleb128 .LVU339
	.uleb128 .LVU398
	.uleb128 .LVU400
.LLST77:
	.byte	0x4
	.uleb128 .LVL86-.Ltext0
	.uleb128 .LVL88-.Ltext0
	.uleb128 0x1
	.byte	0x5a
	.byte	0x4
	.uleb128 .LVL102-.Ltext0
	.uleb128 .LVL104-.Ltext0
	.uleb128 0x1
	.byte	0x5a
	.byte	0
.LVUS79:
	.uleb128 .LVU342
	.uleb128 .LVU345
.LLST79:
	.byte	0x4
	.uleb128 .LVL89-.Ltext0
	.uleb128 .LVL89-.Ltext0
	.uleb128 0x4
	.byte	0xa3
	.uleb128 0x1
	.byte	0x52
	.byte	0x9f
	.byte	0
.LVUS80:
	.uleb128 .LVU342
	.uleb128 .LVU345
.LLST80:
	.byte	0x4
	.uleb128 .LVL89-.Ltext0
	.uleb128 .LVL89-.Ltext0
	.uleb128 0x1
	.byte	0x51
	.byte	0
.LVUS81:
	.uleb128 .LVU344
	.uleb128 .LVU345
.LLST81:
	.byte	0x4
	.uleb128 .LVL89-.Ltext0
	.uleb128 .LVL89-.Ltext0
	.uleb128 0xe
	.byte	0x71
	.sleb128 8
	.byte	0x6
	.byte	0x75
	.sleb128 0
	.byte	0x22
	.byte	0xa6
	.byte	0x8
	.uleb128 0x2e
	.byte	0xa5
	.uleb128 0x12
	.uleb128 0x2e
	.byte	0x1c
	.byte	0x9f
	.byte	0
.LVUS82:
	.uleb128 .LVU353
	.uleb128 .LVU361
.LLST82:
	.byte	0x4
	.uleb128 .LVL91-.Ltext0
	.uleb128 .LVL93-.Ltext0
	.uleb128 0x1
	.byte	0x50
	.byte	0
.LVUS83:
	.uleb128 .LVU353
	.uleb128 .LVU361
.LLST83:
	.byte	0x4
	.uleb128 .LVL91-.Ltext0
	.uleb128 .LVL93-.Ltext0
	.uleb128 0x1
	.byte	0x52
	.byte	0
.LVUS84:
	.uleb128 .LVU357
	.uleb128 .LVU361
.LLST84:
	.byte	0x4
	.uleb128 .LVL92-.Ltext0
	.uleb128 .LVL93-.Ltext0
	.uleb128 0x1
	.byte	0x5e
	.byte	0
.LVUS86:
	.uleb128 .LVU370
	.uleb128 .LVU378
.LLST86:
	.byte	0x4
	.uleb128 .LVL96-.Ltext0
	.uleb128 .LVL98-.Ltext0
	.uleb128 0x1
	.byte	0x53
	.byte	0
.LVUS87:
	.uleb128 .LVU370
	.uleb128 .LVU378
.LLST87:
	.byte	0x4
	.uleb128 .LVL96-.Ltext0
	.uleb128 .LVL98-.Ltext0
	.uleb128 0x1
	.byte	0x51
	.byte	0
.LVUS88:
	.uleb128 .LVU374
	.uleb128 .LVU378
.LLST88:
	.byte	0x4
	.uleb128 .LVL97-.Ltext0
	.uleb128 .LVL98-.Ltext0
	.uleb128 0x1
	.byte	0x5e
	.byte	0
.LVUS89:
	.uleb128 .LVU405
	.uleb128 .LVU433
.LLST89:
	.byte	0x4
	.uleb128 .LVL106-.Ltext0
	.uleb128 .LVL116-.Ltext0
	.uleb128 0x1
	.byte	0x52
	.byte	0
.LVUS90:
	.uleb128 .LVU405
	.uleb128 .LVU409
	.uleb128 .LVU409
	.uleb128 .LVU420
	.uleb128 .LVU420
	.uleb128 .LVU433
.LLST90:
	.byte	0x4
	.uleb128 .LVL106-.Ltext0
	.uleb128 .LVL108-.Ltext0
	.uleb128 0x1
	.byte	0x51
	.byte	0x4
	.uleb128 .LVL108-.Ltext0
	.uleb128 .LVL112-.Ltext0
	.uleb128 0x1
	.byte	0x58
	.byte	0x4
	.uleb128 .LVL112-.Ltext0
	.uleb128 .LVL116-.Ltext0
	.uleb128 0x1
	.byte	0x5b
	.byte	0
.LVUS91:
	.uleb128 .LVU405
	.uleb128 .LVU408
	.uleb128 .LVU408
	.uleb128 .LVU409
	.uleb128 .LVU409
	.uleb128 .LVU420
	.uleb128 .LVU420
	.uleb128 .LVU433
.LLST91:
	.byte	0x4
	.uleb128 .LVL106-.Ltext0
	.uleb128 .LVL107-.Ltext0
	.uleb128 0x1
	.byte	0x50
	.byte	0x4
	.uleb128 .LVL107-.Ltext0
	.uleb128 .LVL108-.Ltext0
	.uleb128 0x3
	.byte	0x71
	.sleb128 -1
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL108-.Ltext0
	.uleb128 .LVL112-.Ltext0
	.uleb128 0x3
	.byte	0x78
	.sleb128 -1
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL112-.Ltext0
	.uleb128 .LVL116-.Ltext0
	.uleb128 0x3
	.byte	0x7b
	.sleb128 -1
	.byte	0x9f
	.byte	0
.LVUS92:
	.uleb128 .LVU405
	.uleb128 .LVU415
	.uleb128 .LVU415
	.uleb128 .LVU421
	.uleb128 .LVU421
	.uleb128 .LVU433
.LLST92:
	.byte	0x4
	.uleb128 .LVL106-.Ltext0
	.uleb128 .LVL110-.Ltext0
	.uleb128 0x1
	.byte	0x55
	.byte	0x4
	.uleb128 .LVL110-.Ltext0
	.uleb128 .LVL113-.Ltext0
	.uleb128 0x1
	.byte	0x59
	.byte	0x4
	.uleb128 .LVL113-.Ltext0
	.uleb128 .LVL116-.Ltext0
	.uleb128 0x4
	.byte	0xa3
	.uleb128 0x1
	.byte	0x55
	.byte	0x9f
	.byte	0
.LVUS94:
	.uleb128 .LVU416
	.uleb128 .LVU419
.LLST94:
	.byte	0x4
	.uleb128 .LVL111-.Ltext0
	.uleb128 .LVL111-.Ltext0
	.uleb128 0x1
	.byte	0x52
	.byte	0
.LVUS95:
	.uleb128 .LVU416
	.uleb128 .LVU419
.LLST95:
	.byte	0x4
	.uleb128 .LVL111-.Ltext0
	.uleb128 .LVL111-.Ltext0
	.uleb128 0x1
	.byte	0x55
	.byte	0
.LVUS96:
	.uleb128 .LVU416
	.uleb128 .LVU419
.LLST96:
	.byte	0x4
	.uleb128 .LVL111-.Ltext0
	.uleb128 .LVL111-.Ltext0
	.uleb128 0x1
	.byte	0x50
	.byte	0
.LVUS97:
	.uleb128 .LVU418
	.uleb128 .LVU419
.LLST97:
	.byte	0x4
	.uleb128 .LVL111-.Ltext0
	.uleb128 .LVL111-.Ltext0
	.uleb128 0x18
	.byte	0x72
	.sleb128 0
	.byte	0x33
	.byte	0x24
	.byte	0x70
	.sleb128 8
	.byte	0x6
	.byte	0x22
	.byte	0xa6
	.byte	0x8
	.uleb128 0x2e
	.byte	0x72
	.sleb128 0
	.byte	0x33
	.byte	0x24
	.byte	0x75
	.sleb128 8
	.byte	0x6
	.byte	0x22
	.byte	0xa6
	.byte	0x8
	.uleb128 0x2e
	.byte	0x1c
	.byte	0x9f
	.byte	0
.LVUS98:
	.uleb128 .LVU424
	.uleb128 .LVU431
.LLST98:
	.byte	0x4
	.uleb128 .LVL114-.Ltext0
	.uleb128 .LVL115-.Ltext0
	.uleb128 0x1
	.byte	0x51
	.byte	0
.LVUS99:
	.uleb128 .LVU424
	.uleb128 .LVU431
.LLST99:
	.byte	0x4
	.uleb128 .LVL114-.Ltext0
	.uleb128 .LVL115-.Ltext0
	.uleb128 0x1
	.byte	0x54
	.byte	0
.LVUS100:
	.uleb128 .LVU427
	.uleb128 .LVU431
.LLST100:
	.byte	0x4
	.uleb128 .LVL114-.Ltext0
	.uleb128 .LVL115-.Ltext0
	.uleb128 0x1
	.byte	0x50
	.byte	0
.LVUS101:
	.uleb128 0
	.uleb128 .LVU454
	.uleb128 .LVU454
	.uleb128 .LVU475
	.uleb128 .LVU475
	.uleb128 .LVU478
	.uleb128 .LVU478
	.uleb128 0
.LLST101:
	.byte	0x4
	.uleb128 .LVL118-.Ltext0
	.uleb128 .LVL124-1-.Ltext0
	.uleb128 0x1
	.byte	0x55
	.byte	0x4
	.uleb128 .LVL124-1-.Ltext0
	.uleb128 .LVL130-.Ltext0
	.uleb128 0x3
	.byte	0x91
	.sleb128 -80
	.byte	0x4
	.uleb128 .LVL130-.Ltext0
	.uleb128 .LVL133-.Ltext0
	.uleb128 0x4
	.byte	0xa3
	.uleb128 0x1
	.byte	0x55
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL133-.Ltext0
	.uleb128 .LFE63-.Ltext0
	.uleb128 0x1
	.byte	0x55
	.byte	0
.LVUS102:
	.uleb128 0
	.uleb128 .LVU454
	.uleb128 .LVU454
	.uleb128 .LVU475
	.uleb128 .LVU475
	.uleb128 .LVU478
	.uleb128 .LVU478
	.uleb128 0
.LLST102:
	.byte	0x4
	.uleb128 .LVL118-.Ltext0
	.uleb128 .LVL124-1-.Ltext0
	.uleb128 0x1
	.byte	0x54
	.byte	0x4
	.uleb128 .LVL124-1-.Ltext0
	.uleb128 .LVL130-.Ltext0
	.uleb128 0x3
	.byte	0x91
	.sleb128 -72
	.byte	0x4
	.uleb128 .LVL130-.Ltext0
	.uleb128 .LVL133-.Ltext0
	.uleb128 0x4
	.byte	0xa3
	.uleb128 0x1
	.byte	0x54
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL133-.Ltext0
	.uleb128 .LFE63-.Ltext0
	.uleb128 0x1
	.byte	0x54
	.byte	0
.LVUS103:
	.uleb128 0
	.uleb128 .LVU443
	.uleb128 .LVU443
	.uleb128 .LVU476
	.uleb128 .LVU476
	.uleb128 0
.LLST103:
	.byte	0x4
	.uleb128 .LVL118-.Ltext0
	.uleb128 .LVL120-.Ltext0
	.uleb128 0x1
	.byte	0x51
	.byte	0x4
	.uleb128 .LVL120-.Ltext0
	.uleb128 .LVL131-.Ltext0
	.uleb128 0x1
	.byte	0x5f
	.byte	0x4
	.uleb128 .LVL131-.Ltext0
	.uleb128 .LFE63-.Ltext0
	.uleb128 0x4
	.byte	0xa3
	.uleb128 0x1
	.byte	0x51
	.byte	0x9f
	.byte	0
.LVUS104:
	.uleb128 0
	.uleb128 .LVU452
	.uleb128 .LVU452
	.uleb128 .LVU475
	.uleb128 .LVU475
	.uleb128 .LVU478
	.uleb128 .LVU478
	.uleb128 0
.LLST104:
	.byte	0x4
	.uleb128 .LVL118-.Ltext0
	.uleb128 .LVL122-.Ltext0
	.uleb128 0x1
	.byte	0x52
	.byte	0x4
	.uleb128 .LVL122-.Ltext0
	.uleb128 .LVL130-.Ltext0
	.uleb128 0x1
	.byte	0x5d
	.byte	0x4
	.uleb128 .LVL130-.Ltext0
	.uleb128 .LVL133-.Ltext0
	.uleb128 0x4
	.byte	0xa3
	.uleb128 0x1
	.byte	0x52
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL133-.Ltext0
	.uleb128 .LFE63-.Ltext0
	.uleb128 0x1
	.byte	0x52
	.byte	0
.LVUS105:
	.uleb128 0
	.uleb128 .LVU454
	.uleb128 .LVU454
	.uleb128 .LVU475
	.uleb128 .LVU475
	.uleb128 .LVU478
	.uleb128 .LVU478
	.uleb128 0
.LLST105:
	.byte	0x4
	.uleb128 .LVL118-.Ltext0
	.uleb128 .LVL124-1-.Ltext0
	.uleb128 0x1
	.byte	0x58
	.byte	0x4
	.uleb128 .LVL124-1-.Ltext0
	.uleb128 .LVL130-.Ltext0
	.uleb128 0x1
	.byte	0x56
	.byte	0x4
	.uleb128 .LVL130-.Ltext0
	.uleb128 .LVL133-.Ltext0
	.uleb128 0x4
	.byte	0xa3
	.uleb128 0x1
	.byte	0x58
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL133-.Ltext0
	.uleb128 .LFE63-.Ltext0
	.uleb128 0x1
	.byte	0x56
	.byte	0
.LVUS106:
	.uleb128 .LVU438
	.uleb128 .LVU462
	.uleb128 .LVU462
	.uleb128 .LVU475
	.uleb128 .LVU477
	.uleb128 .LVU480
	.uleb128 .LVU480
	.uleb128 0
.LLST106:
	.byte	0x4
	.uleb128 .LVL119-.Ltext0
	.uleb128 .LVL127-.Ltext0
	.uleb128 0x2
	.byte	0x30
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL127-.Ltext0
	.uleb128 .LVL130-.Ltext0
	.uleb128 0x1
	.byte	0x53
	.byte	0x4
	.uleb128 .LVL132-.Ltext0
	.uleb128 .LVL134-.Ltext0
	.uleb128 0x2
	.byte	0x30
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL134-.Ltext0
	.uleb128 .LFE63-.Ltext0
	.uleb128 0x1
	.byte	0x53
	.byte	0
.LVUS107:
	.uleb128 .LVU445
	.uleb128 .LVU453
	.uleb128 .LVU453
	.uleb128 .LVU454
	.uleb128 .LVU454
	.uleb128 .LVU475
	.uleb128 .LVU478
	.uleb128 0
.LLST107:
	.byte	0x4
	.uleb128 .LVL121-.Ltext0
	.uleb128 .LVL123-.Ltext0
	.uleb128 0x1
	.byte	0x51
	.byte	0x4
	.uleb128 .LVL123-.Ltext0
	.uleb128 .LVL124-1-.Ltext0
	.uleb128 0x1
	.byte	0x52
	.byte	0x4
	.uleb128 .LVL124-1-.Ltext0
	.uleb128 .LVL130-.Ltext0
	.uleb128 0x1
	.byte	0x5e
	.byte	0x4
	.uleb128 .LVL133-.Ltext0
	.uleb128 .LFE63-.Ltext0
	.uleb128 0x1
	.byte	0x51
	.byte	0
.LVUS108:
	.uleb128 .LVU446
	.uleb128 .LVU455
	.uleb128 .LVU455
	.uleb128 .LVU460
	.uleb128 .LVU460
	.uleb128 .LVU475
	.uleb128 .LVU478
	.uleb128 0
.LLST108:
	.byte	0x4
	.uleb128 .LVL121-.Ltext0
	.uleb128 .LVL125-.Ltext0
	.uleb128 0x3
	.byte	0x9
	.byte	0xff
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL125-.Ltext0
	.uleb128 .LVL126-.Ltext0
	.uleb128 0x1
	.byte	0x50
	.byte	0x4
	.uleb128 .LVL126-.Ltext0
	.uleb128 .LVL130-.Ltext0
	.uleb128 0x1
	.byte	0x5c
	.byte	0x4
	.uleb128 .LVL133-.Ltext0
	.uleb128 .LFE63-.Ltext0
	.uleb128 0x3
	.byte	0x9
	.byte	0xff
	.byte	0x9f
	.byte	0
.LVUS109:
	.uleb128 0
	.uleb128 .LVU511
	.uleb128 .LVU511
	.uleb128 .LVU531
	.uleb128 .LVU531
	.uleb128 .LVU534
	.uleb128 .LVU534
	.uleb128 0
.LLST109:
	.byte	0x4
	.uleb128 .LVL135-.Ltext0
	.uleb128 .LVL140-1-.Ltext0
	.uleb128 0x1
	.byte	0x55
	.byte	0x4
	.uleb128 .LVL140-1-.Ltext0
	.uleb128 .LVL146-.Ltext0
	.uleb128 0x1
	.byte	0x5c
	.byte	0x4
	.uleb128 .LVL146-.Ltext0
	.uleb128 .LVL149-.Ltext0
	.uleb128 0x4
	.byte	0xa3
	.uleb128 0x1
	.byte	0x55
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL149-.Ltext0
	.uleb128 .LFE71-.Ltext0
	.uleb128 0x1
	.byte	0x55
	.byte	0
.LVUS110:
	.uleb128 0
	.uleb128 .LVU511
	.uleb128 .LVU511
	.uleb128 .LVU531
	.uleb128 .LVU531
	.uleb128 .LVU534
	.uleb128 .LVU534
	.uleb128 0
.LLST110:
	.byte	0x4
	.uleb128 .LVL135-.Ltext0
	.uleb128 .LVL140-1-.Ltext0
	.uleb128 0x1
	.byte	0x54
	.byte	0x4
	.uleb128 .LVL140-1-.Ltext0
	.uleb128 .LVL146-.Ltext0
	.uleb128 0x1
	.byte	0x5d
	.byte	0x4
	.uleb128 .LVL146-.Ltext0
	.uleb128 .LVL149-.Ltext0
	.uleb128 0x4
	.byte	0xa3
	.uleb128 0x1
	.byte	0x54
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL149-.Ltext0
	.uleb128 .LFE71-.Ltext0
	.uleb128 0x1
	.byte	0x54
	.byte	0
.LVUS111:
	.uleb128 0
	.uleb128 .LVU499
	.uleb128 .LVU499
	.uleb128 .LVU532
	.uleb128 .LVU532
	.uleb128 0
.LLST111:
	.byte	0x4
	.uleb128 .LVL135-.Ltext0
	.uleb128 .LVL136-.Ltext0
	.uleb128 0x1
	.byte	0x51
	.byte	0x4
	.uleb128 .LVL136-.Ltext0
	.uleb128 .LVL147-.Ltext0
	.uleb128 0x1
	.byte	0x56
	.byte	0x4
	.uleb128 .LVL147-.Ltext0
	.uleb128 .LFE71-.Ltext0
	.uleb128 0x4
	.byte	0xa3
	.uleb128 0x1
	.byte	0x51
	.byte	0x9f
	.byte	0
.LVUS112:
	.uleb128 0
	.uleb128 .LVU507
	.uleb128 .LVU507
	.uleb128 .LVU531
	.uleb128 .LVU531
	.uleb128 .LVU534
	.uleb128 .LVU534
	.uleb128 0
.LLST112:
	.byte	0x4
	.uleb128 .LVL135-.Ltext0
	.uleb128 .LVL138-.Ltext0
	.uleb128 0x1
	.byte	0x52
	.byte	0x4
	.uleb128 .LVL138-.Ltext0
	.uleb128 .LVL146-.Ltext0
	.uleb128 0x1
	.byte	0x5e
	.byte	0x4
	.uleb128 .LVL146-.Ltext0
	.uleb128 .LVL149-.Ltext0
	.uleb128 0x4
	.byte	0xa3
	.uleb128 0x1
	.byte	0x52
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL149-.Ltext0
	.uleb128 .LFE71-.Ltext0
	.uleb128 0x1
	.byte	0x52
	.byte	0
.LVUS113:
	.uleb128 .LVU494
	.uleb128 .LVU518
	.uleb128 .LVU518
	.uleb128 .LVU531
	.uleb128 .LVU533
	.uleb128 .LVU536
	.uleb128 .LVU536
	.uleb128 0
.LLST113:
	.byte	0x4
	.uleb128 .LVL135-.Ltext0
	.uleb128 .LVL143-.Ltext0
	.uleb128 0x2
	.byte	0x30
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL143-.Ltext0
	.uleb128 .LVL146-.Ltext0
	.uleb128 0x1
	.byte	0x53
	.byte	0x4
	.uleb128 .LVL148-.Ltext0
	.uleb128 .LVL150-.Ltext0
	.uleb128 0x2
	.byte	0x30
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL150-.Ltext0
	.uleb128 .LFE71-.Ltext0
	.uleb128 0x1
	.byte	0x53
	.byte	0
.LVUS114:
	.uleb128 .LVU501
	.uleb128 .LVU508
	.uleb128 .LVU508
	.uleb128 .LVU511
	.uleb128 .LVU511
	.uleb128 .LVU531
	.uleb128 .LVU534
	.uleb128 0
.LLST114:
	.byte	0x4
	.uleb128 .LVL137-.Ltext0
	.uleb128 .LVL139-.Ltext0
	.uleb128 0x1
	.byte	0x51
	.byte	0x4
	.uleb128 .LVL139-.Ltext0
	.uleb128 .LVL140-1-.Ltext0
	.uleb128 0x1
	.byte	0x52
	.byte	0x4
	.uleb128 .LVL140-1-.Ltext0
	.uleb128 .LVL146-.Ltext0
	.uleb128 0x1
	.byte	0x5f
	.byte	0x4
	.uleb128 .LVL149-.Ltext0
	.uleb128 .LFE71-.Ltext0
	.uleb128 0x1
	.byte	0x51
	.byte	0
.LVUS115:
	.uleb128 .LVU502
	.uleb128 .LVU512
	.uleb128 .LVU512
	.uleb128 .LVU516
	.uleb128 .LVU516
	.uleb128 .LVU520
	.uleb128 .LVU520
	.uleb128 .LVU531
	.uleb128 .LVU534
	.uleb128 0
.LLST115:
	.byte	0x4
	.uleb128 .LVL137-.Ltext0
	.uleb128 .LVL141-.Ltext0
	.uleb128 0x3
	.byte	0x9
	.byte	0xff
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL141-.Ltext0
	.uleb128 .LVL142-.Ltext0
	.uleb128 0x1
	.byte	0x50
	.byte	0x4
	.uleb128 .LVL142-.Ltext0
	.uleb128 .LVL144-1-.Ltext0
	.uleb128 0x1
	.byte	0x59
	.byte	0x4
	.uleb128 .LVL144-1-.Ltext0
	.uleb128 .LVL146-.Ltext0
	.uleb128 0x3
	.byte	0x91
	.sleb128 -68
	.byte	0x4
	.uleb128 .LVL149-.Ltext0
	.uleb128 .LFE71-.Ltext0
	.uleb128 0x3
	.byte	0x9
	.byte	0xff
	.byte	0x9f
	.byte	0
.LVUS116:
	.uleb128 0
	.uleb128 .LVU567
	.uleb128 .LVU567
	.uleb128 .LVU587
	.uleb128 .LVU587
	.uleb128 .LVU590
	.uleb128 .LVU590
	.uleb128 0
.LLST116:
	.byte	0x4
	.uleb128 .LVL151-.Ltext0
	.uleb128 .LVL156-1-.Ltext0
	.uleb128 0x1
	.byte	0x55
	.byte	0x4
	.uleb128 .LVL156-1-.Ltext0
	.uleb128 .LVL162-.Ltext0
	.uleb128 0x1
	.byte	0x5c
	.byte	0x4
	.uleb128 .LVL162-.Ltext0
	.uleb128 .LVL165-.Ltext0
	.uleb128 0x4
	.byte	0xa3
	.uleb128 0x1
	.byte	0x55
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL165-.Ltext0
	.uleb128 .LFE72-.Ltext0
	.uleb128 0x1
	.byte	0x55
	.byte	0
.LVUS117:
	.uleb128 0
	.uleb128 .LVU567
	.uleb128 .LVU567
	.uleb128 .LVU587
	.uleb128 .LVU587
	.uleb128 .LVU590
	.uleb128 .LVU590
	.uleb128 0
.LLST117:
	.byte	0x4
	.uleb128 .LVL151-.Ltext0
	.uleb128 .LVL156-1-.Ltext0
	.uleb128 0x1
	.byte	0x54
	.byte	0x4
	.uleb128 .LVL156-1-.Ltext0
	.uleb128 .LVL162-.Ltext0
	.uleb128 0x1
	.byte	0x5d
	.byte	0x4
	.uleb128 .LVL162-.Ltext0
	.uleb128 .LVL165-.Ltext0
	.uleb128 0x4
	.byte	0xa3
	.uleb128 0x1
	.byte	0x54
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL165-.Ltext0
	.uleb128 .LFE72-.Ltext0
	.uleb128 0x1
	.byte	0x54
	.byte	0
.LVUS118:
	.uleb128 0
	.uleb128 .LVU555
	.uleb128 .LVU555
	.uleb128 .LVU588
	.uleb128 .LVU588
	.uleb128 0
.LLST118:
	.byte	0x4
	.uleb128 .LVL151-.Ltext0
	.uleb128 .LVL152-.Ltext0
	.uleb128 0x1
	.byte	0x51
	.byte	0x4
	.uleb128 .LVL152-.Ltext0
	.uleb128 .LVL163-.Ltext0
	.uleb128 0x1
	.byte	0x56
	.byte	0x4
	.uleb128 .LVL163-.Ltext0
	.uleb128 .LFE72-.Ltext0
	.uleb128 0x4
	.byte	0xa3
	.uleb128 0x1
	.byte	0x51
	.byte	0x9f
	.byte	0
.LVUS119:
	.uleb128 0
	.uleb128 .LVU563
	.uleb128 .LVU563
	.uleb128 .LVU587
	.uleb128 .LVU587
	.uleb128 .LVU590
	.uleb128 .LVU590
	.uleb128 0
.LLST119:
	.byte	0x4
	.uleb128 .LVL151-.Ltext0
	.uleb128 .LVL154-.Ltext0
	.uleb128 0x1
	.byte	0x52
	.byte	0x4
	.uleb128 .LVL154-.Ltext0
	.uleb128 .LVL162-.Ltext0
	.uleb128 0x1
	.byte	0x5e
	.byte	0x4
	.uleb128 .LVL162-.Ltext0
	.uleb128 .LVL165-.Ltext0
	.uleb128 0x4
	.byte	0xa3
	.uleb128 0x1
	.byte	0x52
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL165-.Ltext0
	.uleb128 .LFE72-.Ltext0
	.uleb128 0x1
	.byte	0x52
	.byte	0
.LVUS120:
	.uleb128 .LVU550
	.uleb128 .LVU574
	.uleb128 .LVU574
	.uleb128 .LVU587
	.uleb128 .LVU589
	.uleb128 .LVU592
	.uleb128 .LVU592
	.uleb128 0
.LLST120:
	.byte	0x4
	.uleb128 .LVL151-.Ltext0
	.uleb128 .LVL159-.Ltext0
	.uleb128 0x2
	.byte	0x30
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL159-.Ltext0
	.uleb128 .LVL162-.Ltext0
	.uleb128 0x1
	.byte	0x53
	.byte	0x4
	.uleb128 .LVL164-.Ltext0
	.uleb128 .LVL166-.Ltext0
	.uleb128 0x2
	.byte	0x30
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL166-.Ltext0
	.uleb128 .LFE72-.Ltext0
	.uleb128 0x1
	.byte	0x53
	.byte	0
.LVUS121:
	.uleb128 .LVU557
	.uleb128 .LVU564
	.uleb128 .LVU564
	.uleb128 .LVU567
	.uleb128 .LVU567
	.uleb128 .LVU587
	.uleb128 .LVU590
	.uleb128 0
.LLST121:
	.byte	0x4
	.uleb128 .LVL153-.Ltext0
	.uleb128 .LVL155-.Ltext0
	.uleb128 0x1
	.byte	0x51
	.byte	0x4
	.uleb128 .LVL155-.Ltext0
	.uleb128 .LVL156-1-.Ltext0
	.uleb128 0x1
	.byte	0x52
	.byte	0x4
	.uleb128 .LVL156-1-.Ltext0
	.uleb128 .LVL162-.Ltext0
	.uleb128 0x1
	.byte	0x5f
	.byte	0x4
	.uleb128 .LVL165-.Ltext0
	.uleb128 .LFE72-.Ltext0
	.uleb128 0x1
	.byte	0x51
	.byte	0
.LVUS122:
	.uleb128 .LVU558
	.uleb128 .LVU568
	.uleb128 .LVU568
	.uleb128 .LVU572
	.uleb128 .LVU572
	.uleb128 .LVU576
	.uleb128 .LVU576
	.uleb128 .LVU587
	.uleb128 .LVU590
	.uleb128 0
.LLST122:
	.byte	0x4
	.uleb128 .LVL153-.Ltext0
	.uleb128 .LVL157-.Ltext0
	.uleb128 0x3
	.byte	0x9
	.byte	0xff
	.byte	0x9f
	.byte	0x4
	.uleb128 .LVL157-.Ltext0
	.uleb128 .LVL158-.Ltext0
	.uleb128 0x1
	.byte	0x50
	.byte	0x4
	.uleb128 .LVL158-.Ltext0
	.uleb128 .LVL160-1-.Ltext0
	.uleb128 0x1
	.byte	0x58
	.byte	0x4
	.uleb128 .LVL160-1-.Ltext0
	.uleb128 .LVL162-.Ltext0
	.uleb128 0x3
	.byte	0x91
	.sleb128 -68
	.byte	0x4
	.uleb128 .LVL165-.Ltext0
	.uleb128 .LFE72-.Ltext0
	.uleb128 0x3
	.byte	0x9
	.byte	0xff
	.byte	0x9f
	.byte	0
.Ldebug_loc3:
	.section	.debug_aranges,"",@progbits
	.long	0x2c
	.value	0x2
	.long	.Ldebug_info0
	.byte	0x8
	.byte	0
	.value	0
	.value	0
	.quad	.Ltext0
	.quad	.Letext0-.Ltext0
	.quad	0
	.quad	0
	.section	.debug_rnglists,"",@progbits
.Ldebug_ranges0:
	.long	.Ldebug_ranges3-.Ldebug_ranges2
.Ldebug_ranges2:
	.value	0x5
	.byte	0x8
	.byte	0
	.long	0
.LLRL11:
	.byte	0x4
	.uleb128 .LBB66-.Ltext0
	.uleb128 .LBE66-.Ltext0
	.byte	0x4
	.uleb128 .LBB70-.Ltext0
	.uleb128 .LBE70-.Ltext0
	.byte	0x4
	.uleb128 .LBB71-.Ltext0
	.uleb128 .LBE71-.Ltext0
	.byte	0x4
	.uleb128 .LBB72-.Ltext0
	.uleb128 .LBE72-.Ltext0
	.byte	0
.LLRL13:
	.byte	0x4
	.uleb128 .LBB67-.Ltext0
	.uleb128 .LBE67-.Ltext0
	.byte	0x4
	.uleb128 .LBB68-.Ltext0
	.uleb128 .LBE68-.Ltext0
	.byte	0x4
	.uleb128 .LBB69-.Ltext0
	.uleb128 .LBE69-.Ltext0
	.byte	0
.LLRL27:
	.byte	0x4
	.uleb128 .LBB80-.Ltext0
	.uleb128 .LBE80-.Ltext0
	.byte	0x4
	.uleb128 .LBB81-.Ltext0
	.uleb128 .LBE81-.Ltext0
	.byte	0x4
	.uleb128 .LBB82-.Ltext0
	.uleb128 .LBE82-.Ltext0
	.byte	0
.LLRL31:
	.byte	0x4
	.uleb128 .LBB83-.Ltext0
	.uleb128 .LBE83-.Ltext0
	.byte	0x4
	.uleb128 .LBB88-.Ltext0
	.uleb128 .LBE88-.Ltext0
	.byte	0x4
	.uleb128 .LBB89-.Ltext0
	.uleb128 .LBE89-.Ltext0
	.byte	0x4
	.uleb128 .LBB90-.Ltext0
	.uleb128 .LBE90-.Ltext0
	.byte	0
.LLRL36:
	.byte	0x4
	.uleb128 .LBB96-.Ltext0
	.uleb128 .LBE96-.Ltext0
	.byte	0x4
	.uleb128 .LBB100-.Ltext0
	.uleb128 .LBE100-.Ltext0
	.byte	0x4
	.uleb128 .LBB101-.Ltext0
	.uleb128 .LBE101-.Ltext0
	.byte	0
.LLRL44:
	.byte	0x4
	.uleb128 .LBB114-.Ltext0
	.uleb128 .LBE114-.Ltext0
	.byte	0x4
	.uleb128 .LBB118-.Ltext0
	.uleb128 .LBE118-.Ltext0
	.byte	0x4
	.uleb128 .LBB119-.Ltext0
	.uleb128 .LBE119-.Ltext0
	.byte	0
.LLRL51:
	.byte	0x4
	.uleb128 .LBB128-.Ltext0
	.uleb128 .LBE128-.Ltext0
	.byte	0x4
	.uleb128 .LBB132-.Ltext0
	.uleb128 .LBE132-.Ltext0
	.byte	0x4
	.uleb128 .LBB133-.Ltext0
	.uleb128 .LBE133-.Ltext0
	.byte	0
.LLRL58:
	.byte	0x4
	.uleb128 .LBB136-.Ltext0
	.uleb128 .LBE136-.Ltext0
	.byte	0x4
	.uleb128 .LBB139-.Ltext0
	.uleb128 .LBE139-.Ltext0
	.byte	0
.LLRL67:
	.byte	0x4
	.uleb128 .LBB157-.Ltext0
	.uleb128 .LBE157-.Ltext0
	.byte	0x4
	.uleb128 .LBB184-.Ltext0
	.uleb128 .LBE184-.Ltext0
	.byte	0x4
	.uleb128 .LBB185-.Ltext0
	.uleb128 .LBE185-.Ltext0
	.byte	0x4
	.uleb128 .LBB186-.Ltext0
	.uleb128 .LBE186-.Ltext0
	.byte	0
.LLRL69:
	.byte	0x4
	.uleb128 .LBB158-.Ltext0
	.uleb128 .LBE158-.Ltext0
	.byte	0x4
	.uleb128 .LBB181-.Ltext0
	.uleb128 .LBE181-.Ltext0
	.byte	0x4
	.uleb128 .LBB182-.Ltext0
	.uleb128 .LBE182-.Ltext0
	.byte	0x4
	.uleb128 .LBB183-.Ltext0
	.uleb128 .LBE183-.Ltext0
	.byte	0
.LLRL76:
	.byte	0x4
	.uleb128 .LBB160-.Ltext0
	.uleb128 .LBE160-.Ltext0
	.byte	0x4
	.uleb128 .LBB171-.Ltext0
	.uleb128 .LBE171-.Ltext0
	.byte	0x4
	.uleb128 .LBB172-.Ltext0
	.uleb128 .LBE172-.Ltext0
	.byte	0x4
	.uleb128 .LBB177-.Ltext0
	.uleb128 .LBE177-.Ltext0
	.byte	0
.LLRL78:
	.byte	0x4
	.uleb128 .LBB161-.Ltext0
	.uleb128 .LBE161-.Ltext0
	.byte	0x4
	.uleb128 .LBB166-.Ltext0
	.uleb128 .LBE166-.Ltext0
	.byte	0x4
	.uleb128 .LBB167-.Ltext0
	.uleb128 .LBE167-.Ltext0
	.byte	0x4
	.uleb128 .LBB168-.Ltext0
	.uleb128 .LBE168-.Ltext0
	.byte	0
.LLRL85:
	.byte	0x4
	.uleb128 .LBB173-.Ltext0
	.uleb128 .LBE173-.Ltext0
	.byte	0x4
	.uleb128 .LBB176-.Ltext0
	.uleb128 .LBE176-.Ltext0
	.byte	0
.LLRL93:
	.byte	0x4
	.uleb128 .LBB189-.Ltext0
	.uleb128 .LBE189-.Ltext0
	.byte	0x4
	.uleb128 .LBB192-.Ltext0
	.uleb128 .LBE192-.Ltext0
	.byte	0
.LLRL131:
	.byte	0x4
	.uleb128 .LBB195-.Ltext0
	.uleb128 .LBE195-.Ltext0
	.byte	0x4
	.uleb128 .LBB210-.Ltext0
	.uleb128 .LBE210-.Ltext0
	.byte	0x4
	.uleb128 .LBB215-.Ltext0
	.uleb128 .LBE215-.Ltext0
	.byte	0x4
	.uleb128 .LBB216-.Ltext0
	.uleb128 .LBE216-.Ltext0
	.byte	0x4
	.uleb128 .LBB217-.Ltext0
	.uleb128 .LBE217-.Ltext0
	.byte	0
.LLRL134:
	.byte	0x4
	.uleb128 .LBB197-.Ltext0
	.uleb128 .LBE197-.Ltext0
	.byte	0x4
	.uleb128 .LBB201-.Ltext0
	.uleb128 .LBE201-.Ltext0
	.byte	0x4
	.uleb128 .LBB202-.Ltext0
	.uleb128 .LBE202-.Ltext0
	.byte	0x4
	.uleb128 .LBB203-.Ltext0
	.uleb128 .LBE203-.Ltext0
	.byte	0x4
	.uleb128 .LBB204-.Ltext0
	.uleb128 .LBE204-.Ltext0
	.byte	0x4
	.uleb128 .LBB205-.Ltext0
	.uleb128 .LBE205-.Ltext0
	.byte	0
.LLRL136:
	.byte	0x4
	.uleb128 .LBB198-.Ltext0
	.uleb128 .LBE198-.Ltext0
	.byte	0x4
	.uleb128 .LBB199-.Ltext0
	.uleb128 .LBE199-.Ltext0
	.byte	0x4
	.uleb128 .LBB200-.Ltext0
	.uleb128 .LBE200-.Ltext0
	.byte	0
.LLRL138:
	.byte	0x4
	.uleb128 .LBB211-.Ltext0
	.uleb128 .LBE211-.Ltext0
	.byte	0x4
	.uleb128 .LBB214-.Ltext0
	.uleb128 .LBE214-.Ltext0
	.byte	0
.LLRL149:
	.byte	0x4
	.uleb128 .LBB220-.Ltext0
	.uleb128 .LBE220-.Ltext0
	.byte	0x4
	.uleb128 .LBB227-.Ltext0
	.uleb128 .LBE227-.Ltext0
	.byte	0x4
	.uleb128 .LBB228-.Ltext0
	.uleb128 .LBE228-.Ltext0
	.byte	0x4
	.uleb128 .LBB229-.Ltext0
	.uleb128 .LBE229-.Ltext0
	.byte	0x4
	.uleb128 .LBB230-.Ltext0
	.uleb128 .LBE230-.Ltext0
	.byte	0x4
	.uleb128 .LBB231-.Ltext0
	.uleb128 .LBE231-.Ltext0
	.byte	0
.Ldebug_ranges3:
	.section	.debug_line,"",@progbits
.Ldebug_line0:
	.section	.debug_str,"MS",@progbits,1
.LASF50:
	.string	"right"
.LASF16:
	.string	"uint64_t"
.LASF64:
	.string	"swapHeapNode"
.LASF59:
	.string	"initializePTRS"
.LASF10:
	.string	"size_t"
.LASF48:
	.string	"medianOfNodes"
.LASF25:
	.string	"split_var"
.LASF75:
	.string	"putchar"
.LASF30:
	.string	"allocateHeap"
.LASF61:
	.string	"node_array"
.LASF38:
	.string	"maxk"
.LASF70:
	.string	"memcpy"
.LASF37:
	.string	"kdtree_root"
.LASF73:
	.string	"GNU C17 13.2.0 -mavx2 -mtune=generic -march=x86-64 -ggdb -O3 -fasynchronous-unwind-tables -fstack-protector-strong -fstack-clash-protection -fcf-protection"
.LASF66:
	.string	"swap"
.LASF33:
	.string	"dimensions"
.LASF39:
	.string	"current_distance"
.LASF8:
	.string	"short int"
.LASF17:
	.string	"value"
.LASF29:
	.string	"initHeap"
.LASF26:
	.string	"parent"
.LASF34:
	.string	"root"
.LASF63:
	.string	"swap_kd_node_ptrs"
.LASF42:
	.string	"max_d"
.LASF15:
	.string	"float"
.LASF13:
	.string	"long long int"
.LASF44:
	.string	"hyper_plane_dist"
.LASF31:
	.string	"insertMaxHeap"
.LASF9:
	.string	"long int"
.LASF67:
	.string	"__dest"
.LASF53:
	.string	"high"
.LASF71:
	.string	"printf"
.LASF28:
	.string	"HeapSort"
.LASF27:
	.string	"data_dims"
.LASF51:
	.string	"pivotIndex"
.LASF23:
	.string	"kd_node"
.LASF5:
	.string	"unsigned char"
.LASF68:
	.string	"__src"
.LASF19:
	.string	"heap_node"
.LASF7:
	.string	"signed char"
.LASF14:
	.string	"long long unsigned int"
.LASF55:
	.string	"KNN_sub_tree_search"
.LASF18:
	.string	"array_idx"
.LASF4:
	.string	"unsigned int"
.LASF74:
	.string	"__stack_chk_fail"
.LASF46:
	.string	"start"
.LASF62:
	.string	"initializeKDnodes"
.LASF32:
	.string	"kd_ptrs"
.LASF12:
	.string	"char"
.LASF49:
	.string	"left"
.LASF6:
	.string	"short unsigned int"
.LASF65:
	.string	"euclidean_distance"
.LASF76:
	.string	"__builtin_putchar"
.LASF11:
	.string	"__uint64_t"
.LASF22:
	.string	"data"
.LASF56:
	.string	"printKDnode"
.LASF72:
	.string	"__fmt"
.LASF58:
	.string	"cmpKDnodes"
.LASF3:
	.string	"long unsigned int"
.LASF36:
	.string	"point"
.LASF2:
	.string	"double"
.LASF45:
	.string	"make_tree"
.LASF47:
	.string	"median_idx"
.LASF21:
	.string	"count"
.LASF57:
	.string	"node"
.LASF24:
	.string	"level"
.LASF35:
	.string	"build_tree"
.LASF60:
	.string	"node_ptr_array"
.LASF52:
	.string	"partition"
.LASF69:
	.string	"__len"
.LASF43:
	.string	"__printf_chk"
.LASF41:
	.string	"side"
.LASF54:
	.string	"pivot"
.LASF40:
	.string	"hp_distance"
.LASF20:
	.string	"Heap"
	.section	.debug_line_str,"MS",@progbits,1
.LASF1:
	.string	"/home/francesco/Desktop/dssc/robavaria/cc"
.LASF0:
	.string	"src/kdtree.c"
	.ident	"GCC: (Ubuntu 13.2.0-4ubuntu3) 13.2.0"
	.section	.note.GNU-stack,"",@progbits
	.section	.note.gnu.property,"a"
	.align 8
	.long	1f - 0f
	.long	4f - 1f
	.long	5
0:
	.string	"GNU"
1:
	.align 8
	.long	0xc0000002
	.long	3f - 2f
2:
	.long	0x3
3:
	.align 8
4:
