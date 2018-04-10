.func TI_S
TI_S:
push {lr}
push {r4-r7}
nop
nop
nop
nop
nop
nop
nop
nop
   mov r11,r3
   mov r14,r1
   mov r12,r2
   mov r3,r0
   mov r0,r3
   lsr r1,r0,#1
   mov r0,r14
   and r0,r1
   lsl r0,#3
   mov r2,r12
   and r1,r2
   eor r0,r1
   mov r4,r0
   mov r0,r3
   lsl r1,r0,#12
   lsr r0,#4
   eor r0,r1
   mov r1,r11
   and r0,r1
   mov r5,r0
   mov r0,r5
   lsr r1,r0,#1
   mov r0,r14
   and r0,r1
   lsl r0,#3
   mov r2,r12
   and r1,r2
   eor r0,r1
   mov r6,r0
   mov r0,r5
   lsl r1,r0,#12
   lsr r0,#4
   eor r0,r1
   mov r1,r11
   and r0,r1
   mov r7,r0
   mov r0,r7
   lsr r1,r0,#1
   mov r0,r14
   and r0,r1
   lsl r0,#3
   mov r2,r12
   and r1,r2
   eor r0,r1
   mov r8,r0
   mov r0,r7
   lsl r1,r0,#12
   lsr r0,#4
   eor r0,r1
   mov r1,r11
   and r0,r1
   mov r9,r0
   mov r0,r9
   lsr r1,r0,#1
   mov r0,r14
   and r0,r1
   lsl r0,#3
   mov r2,r12
   and r1,r2
   eor r0,r1
   mov r10,r0
   mov r0,#0
   mov r1,r5
   eor r0,r1
   mov r1,r7
   eor r0,r1
   mov r1,r9
   eor r0,r1
   mov r1,r6
   mov r2,r7
   and r1,r2
   eor r0,r1
   mov r1,r5
   mov r2,r7
   and r1,r2
   eor r0,r1
   mov r1,r4
   mov r2,r7
   and r1,r2
   eor r0,r1
   mov r1,r5
   mov r2,r8
   and r1,r2
   eor r0,r1
   mov r1,r6
   mov r2,r9
   and r1,r2
   eor r0,r1
   mov r1,r3
   mov r2,r7
   and r1,r2
   eor r0,r1
   mov r1,r5
   mov r2,r9
   and r1,r2
   eor r0,r1
   mov r1,r3
   mov r2,r8
   and r1,r2
   eor r0,r1
   mov r1,r4
   mov r2,r9
   and r1,r2
   eor r0,r1
   mov r1,r5
   mov r2,r10
   and r1,r2
   eor r0,r1
   mov r1,r3
   mov r2,r9
   and r1,r2
   eor r0,r1
   mov r1,r3
   mov r2,r10
   and r1,r2
   eor r0,r1
nop
nop
nop
nop
nop
nop
nop
nop
pop {r4-r7}
pop {pc}
.endfunc
.func TI_P
TI_P:
push {lr}
push {r4-r5}
   lsl r3,r0,#4
   and r3,r2 
   lsr r0,#12
   mov r4,#7
   and r0,r4
   lsl r0,#2
   ldr r0,[r1,r0]
   eor r0,r3
pop {r4-r5}
pop {pc}
.endfunc
